#include "../include/meshseal/meshseal.h"
#include "internal/component_classifier.h"
#include "internal/diagnostics.h"
#include "internal/si_count.h"
#include "internal/timer.h"
#include "stages/weld.h"
#include "stages/degenerate.h"
#include "stages/orient.h"
#include "stages/nm_vertex.h"
#include "stages/nm_edge.h"
#include "stages/holes.h"
#include "stages/shells.h"
#include "stages/intersections.h"
#include "stages/thin.h"
#include "stages/soup.h"
#include "stages/fwn_levelset.h"
#include "stages/sliver.h"
#include "stages/nm_patch_remesh.h"
#include "stages/bridge_loops.h"
#include "stages/collapse_nm.h"
#include "stages/tjunction.h"
#include "stages/alpha_wrap.h"
#include "stages/nm_local_repair.h"
#include "stages/nm_carve_refill.h"
#include "stl_io.h"
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <unordered_map>
#include <manifold/manifold.h>

namespace meshseal {

namespace {

// Weld a mesh by exact float32 bit pattern — merges vertices that are
// coincident at STL output precision (what a slicer sees). Used both to
// evaluate a mesh's true defect state and to produce a clean indexed mesh.
Mesh float32_weld(const Mesh& in) {
    Mesh out;
    std::map<std::array<float, 3>, uint32_t> seen;
    std::vector<uint32_t> remap(in.vertices.size());
    for (uint32_t i = 0; i < in.vertices.size(); ++i) {
        std::array<float, 3> key{ static_cast<float>(in.vertices[i][0]),
                                  static_cast<float>(in.vertices[i][1]),
                                  static_cast<float>(in.vertices[i][2]) };
        auto it = seen.find(key);
        if (it == seen.end()) {
            uint32_t idx = static_cast<uint32_t>(out.vertices.size());
            seen.emplace(key, idx);
            out.vertices.push_back({ static_cast<double>(key[0]),
                                     static_cast<double>(key[1]),
                                     static_cast<double>(key[2]) });
            remap[i] = idx;
        } else {
            remap[i] = it->second;
        }
    }
    out.faces.reserve(in.faces.size());
    for (const auto& f : in.faces) {
        uint32_t a = remap[f[0]], b = remap[f[1]], c = remap[f[2]];
        if (a != b && b != c && a != c) out.faces.push_back({ a, b, c });
    }
    return out;
}

} // anonymous namespace

RepairResult repair(const Mesh& mesh, const RepairOptions& opts) {
    RepairResult result;
    result.mesh = mesh;

    // --- Input sanity validation ---
    // Reject corrupt input fast instead of grinding on impossible geometry.
    // A malformed STL (NaN/Inf coordinates, absurd magnitudes, or a mesh
    // that is essentially all degenerate faces) would otherwise spin the
    // pipeline for minutes. Bail immediately with an empty result + note.
    {
        bool bad_coord = false;
        double lo[3] = { 0, 0, 0 }, hi[3] = { 0, 0, 0 };
        bool first = true;
        for (const auto& v : mesh.vertices) {
            for (int k = 0; k < 3; ++k) {
                const double c = v[k];
                if (std::isnan(c) || std::isinf(c)) { bad_coord = true; break; }
                if (first) { lo[k] = hi[k] = c; }
                else { if (c < lo[k]) lo[k] = c; if (c > hi[k]) hi[k] = c; }
            }
            if (bad_coord) break;
            first = false;
        }
        double bbd = 0.0;
        if (!first) {
            const double dx = hi[0]-lo[0], dy = hi[1]-lo[1], dz = hi[2]-lo[2];
            bbd = std::sqrt(dx*dx + dy*dy + dz*dz);
        }
        // A real-world printable mesh is at most a few metres; coordinates
        // are conventionally mm. 1e9 is ~1000 km — anything past that is
        // certainly garbage, not a unit-scale mistake.
        const bool absurd_extent = bbd > 1e9;
        if (bad_coord || absurd_extent) {
            result.notes.push_back(
                bad_coord
                  ? "input rejected: NaN/Inf vertex coordinates (corrupt file)"
                  : "input rejected: implausible geometry extent (corrupt file)");
            result.partial_failure = true;
            result.mesh.vertices.clear();
            result.mesh.faces.clear();
            result.self_intersections = 0;
            return result;
        }
    }

    // "Do no harm" baseline: record whether the input is already a clean
    // watertight solid at slicer (float32) precision. A repair pipeline must
    // never hand back a mesh worse than what it was given — if every stage
    // below somehow ends with a non-watertight result on an input that was
    // already watertight, the final guard restores this.
    const Mesh input_welded = float32_weld(mesh);
    bool input_was_clean = false;
    {
        auto d = internal::compute_diagnostics(input_welded);
        if (d.open_boundary_edges == 0 && d.non_manifold_edges == 0 &&
            !input_welded.vertices.empty()) {
            double lo[3] = { input_welded.vertices[0][0],
                             input_welded.vertices[0][1],
                             input_welded.vertices[0][2] };
            double hi[3] = { lo[0], lo[1], lo[2] };
            for (const auto& v : input_welded.vertices)
                for (int k = 0; k < 3; ++k) {
                    if (v[k] < lo[k]) lo[k] = v[k];
                    if (v[k] > hi[k]) hi[k] = v[k];
                }
            const double dx = hi[0]-lo[0], dy = hi[1]-lo[1], dz = hi[2]-lo[2];
            const double diag_cubed = std::pow(dx*dx+dy*dy+dz*dz, 1.5);
            input_was_clean = std::abs(d.signed_volume) > diag_cubed * 1e-6;
        }
    }

    // --- Optional checkpoint dumps (env var MESHSEAL_DUMP_DIR) ---
    // When set, the repair pipeline writes binary-STL snapshots of the
    // working mesh after each major stage. Files are named
    // `dump_<NN>_<label>.stl`. Intended for offline SI-source analysis:
    // running `si_count.py` against each snapshot identifies which stage
    // introduces self-intersections. Zero overhead when env unset.
    const char* dump_env = std::getenv("MESHSEAL_DUMP_DIR");
    const std::filesystem::path dump_dir = dump_env ? std::filesystem::path(dump_env) : std::filesystem::path();
    int dump_seq = 0;
    auto dump = [&](const char* label, const Mesh& m) {
        if (dump_env == nullptr) return;
        try {
            std::filesystem::create_directories(dump_dir);
            char prefix[16];
            std::snprintf(prefix, sizeof(prefix), "dump_%02d_", dump_seq++);
            write_stl(m, dump_dir / (std::string(prefix) + label + ".stl"));
        } catch (...) {
            // Best-effort diagnostic, never fail the pipeline.
        }
    };
    dump("input", result.mesh);

    // --- profiling: env MESHSEAL_PROFILE → per-stage stderr wall-clock ---
    const bool prof_env = std::getenv("MESHSEAL_PROFILE") != nullptr;
    auto prof_t0 = std::chrono::steady_clock::now();
    auto prof_last = prof_t0;
    auto prof_lap = [&](const char* label) {
        if (!prof_env) return;
        auto now = std::chrono::steady_clock::now();
        double since_last = std::chrono::duration<double>(now - prof_last).count();
        double total = std::chrono::duration<double>(now - prof_t0).count();
        std::fprintf(stderr, "[prof] %-26s %8.2fs  (total %8.2fs)\n",
                     label, since_last, total);
        prof_last = now;
    };

    // --- Stage: weld ---
    if (opts.weld && !result.mesh.vertices.empty()) {
        internal::ScopedTimer t("weld", result.stage_times_ms);
        auto wr = stages::weld(result.mesh, opts.weld_tolerance);
        result.mesh = std::move(wr.mesh);
        result.stages_applied.push_back("weld");
        result.notes.push_back(
            "weld: " + std::to_string(wr.vertices_before) +
            " -> " + std::to_string(wr.vertices_after) + " vertices");
    }

    // --- Stage: degenerate ---
    if (opts.degenerate && !result.mesh.faces.empty()) {
        internal::ScopedTimer t("degenerate", result.stage_times_ms);
        auto dr = stages::remove_degenerate(result.mesh);
        result.mesh = std::move(dr.mesh);
        result.stages_applied.push_back("degenerate");
        if (dr.degenerate_removed > 0)
            result.notes.push_back("degenerate: removed " + std::to_string(dr.degenerate_removed) + " zero-area faces");
        if (dr.duplicate_removed > 0)
            result.notes.push_back("degenerate: removed " + std::to_string(dr.duplicate_removed) + " duplicate faces");
        if (dr.isolated_vertices_removed > 0)
            result.notes.push_back("degenerate: removed " + std::to_string(dr.isolated_vertices_removed) + " isolated vertices");
    }

    // --- Stage: orient ---
    if (opts.orient && !result.mesh.faces.empty()) {
        internal::ScopedTimer t("orient", result.stage_times_ms);
        auto or_ = stages::orient_mesh(result.mesh);
        result.mesh = std::move(or_.mesh);
        result.stages_applied.push_back("orient");
        if (or_.faces_flipped > 0)
            result.notes.push_back("orient: flipped " + std::to_string(or_.faces_flipped) + " faces");
        if (or_.components_flipped > 0)
            result.notes.push_back("orient: " + std::to_string(or_.components_flipped) + " components had winding reversed");
        if (!or_.was_orientable)
            result.notes.push_back("orient: mesh has non-orientable regions (skipped BFS propagation)");
    }

    // --- Stage: nm_vertex ---
    if (opts.non_manifold && !result.mesh.faces.empty()) {
        internal::ScopedTimer t("nm_vertex", result.stage_times_ms);
        auto nmv = stages::split_non_manifold_vertices(result.mesh);
        result.mesh = std::move(nmv.mesh);
        result.stages_applied.push_back("nm_vertex");
        if (nmv.vertices_split > 0)
            result.notes.push_back("nm_vertex: split " + std::to_string(nmv.vertices_split) +
                                   " non-manifold vertices (" + std::to_string(nmv.vertices_added) + " vertices added)");
    }

    // --- Stage: nm_edge ---
    if (opts.non_manifold && !result.mesh.faces.empty()) {
        internal::ScopedTimer t("nm_edge", result.stage_times_ms);
        auto nme = stages::fix_non_manifold_edges(result.mesh);
        result.mesh = std::move(nme.mesh);
        if (nme.had_non_manifold_edges) {
            result.stages_applied.push_back("nm_edge");
            if (nme.faces_removed > 0)
                result.notes.push_back("nm_edge: removed " + std::to_string(nme.faces_removed) +
                                       " fin faces, fixed " + std::to_string(nme.edges_fixed) + " edges");
        }
    }

    // --- Stage: orient (BFS-only re-run after nm_edge) ---
    // Per `algorithm_revision.md` §5a: orient should propagate winding
    // through seams that were NM in the input but became manifold after
    // nm_edge. Re-running the FULL orient_mesh after nm_edge was tried
    // 2026-05-13 and consistently regressed `t_junction_cube` — the
    // signed-volume sign-flip step inverted small components whose signed
    // volume happened to change sign due to face removal, producing an
    // inside-out output.
    //
    // The fix (2026-05-13): orient_mesh now accepts a flag that skips the
    // signed-volume sign-flip step. We re-run only the BFS-winding-
    // propagation pass, which is what the design actually wanted.
    if (opts.orient && !result.mesh.faces.empty()) {
        internal::ScopedTimer t("orient", result.stage_times_ms);
        auto or2 = stages::orient_mesh(result.mesh, /*do_signed_volume_flip=*/false);
        result.mesh = std::move(or2.mesh);
        if (or2.faces_flipped > 0)
            result.notes.push_back("orient (post-nm_edge, BFS-only): flipped " +
                std::to_string(or2.faces_flipped) + " faces");
    }

    dump("pre_classify", result.mesh);

    // =================================================================
    // Phases 5R-10R: Component-Aware Repair (Revised Pipeline Architecture)
    // =================================================================

    // Helper: push a stage name exactly once to stages_applied.
    auto add_stage = [&](const std::string& name) {
        if (std::find(result.stages_applied.begin(), result.stages_applied.end(), name)
                == result.stages_applied.end())
            result.stages_applied.push_back(name);
    };

    // Early-exit for empty mesh.
    if (result.mesh.faces.empty()) {
        auto diag = internal::compute_diagnostics(result.mesh);
        result.component_count = diag.component_count;
        result.watertight = (diag.open_boundary_edges == 0 && diag.non_manifold_edges == 0);
        result.is_volume  = result.watertight && (diag.signed_volume > 0.0);
        return result;
    }

    // --- Cross-shell boundary-loop bridging ---
    // A topological opening shared by two shells — e.g. a head shell with a
    // mouth hole and a separate mouth-cavity shell with its own opening
    // (captain_toad) — appears as two distinct boundary loops in two distinct
    // components. The per-component hole-fill that follows would cap each loop
    // independently, sealing the shared cavity into a flat membrane. Bridge
    // such cross-component loop pairs with a triangle strip FIRST so the
    // cavity stays open. No-op on single-component meshes / unpaired loops.
    if (opts.holes && !result.mesh.faces.empty()) {
        internal::ScopedTimer t("bridge_loops", result.stage_times_ms);
        auto br = stages::bridge_paired_loops(result.mesh);
        if (br.bridges_made > 0) {
            result.mesh = std::move(br.mesh);
            add_stage("bridge_loops");
            result.notes.push_back("bridge_loops: zipped " +
                std::to_string(br.bridges_made) + " cross-shell loop pair(s), added " +
                std::to_string(br.faces_added) + " faces");
            // The zipper strip merges two components into one; re-propagate
            // winding (BFS-only) so the new strip faces are consistent.
            if (opts.orient) {
                auto bor = stages::orient_mesh(result.mesh,
                                               /*do_signed_volume_flip=*/false);
                result.mesh = std::move(bor.mesh);
            }
        }
    }

    // --- Phase 5R: Classify connected components ---
    auto components = internal::classify_components(
        result.mesh,
        opts.soup_planarity_threshold,
        opts.soup_open_ratio_threshold);

    // Post-classification routing (Stage 2 of the design): if the input
    // consists of multiple disjoint ultra-planar OPEN components and
    // nothing else, the whole thing is almost certainly a fragmented
    // volumetric solid — e.g. cube_soup_intersecting's 6 face-pairs.
    // Hole-filling each component independently produces overlapping thin
    // shells that the intersections stage can't reconcile. Demote them
    // all to SOUP so they go through the voxel-oracle + T-T arrangement
    // path together, which can reconstruct the implicit volume.
    //
    // Guard: do NOT demote when even one non-planar component exists —
    // that indicates a real 3-D shape coexisting with the planar pieces
    // (a flat plate sitting next to a cube), in which case the planar
    // pieces are usually meaningful surfaces in their own right, not
    // cube fragments.
    {
        int planar_open_count = 0;
        int non_planar_count = 0;
        int other_class_count = 0;
        for (const auto& ci : components) {
            if (ci.cls == internal::ComponentClass::OPEN) {
                if (ci.planarity_ratio < opts.soup_planarity_threshold) {
                    ++planar_open_count;
                } else {
                    ++non_planar_count;
                }
            } else if (ci.cls == internal::ComponentClass::NO_BOUNDARY) {
                ++other_class_count;
            }
            // SOUP-classified components are fine either way.
        }
        if (planar_open_count >= 2 && non_planar_count == 0 && other_class_count == 0) {
            for (auto& ci : components) {
                if (ci.cls == internal::ComponentClass::OPEN) {
                    ci.cls = internal::ComponentClass::SOUP;
                }
            }
        }
    }

    std::vector<Mesh> closed_queue;  // CLOSED queue
    // Parallel "preserve original on Phase-8R-rescue-collapse" flag. True for
    // components that came from REAL input (NO_BOUNDARY components, OPEN→close
    // results) — geometry the user wants kept. False for synthetic components
    // (FWN+LevelSet output, soup-aggregation reconstructions) — those have no
    // user-meaningful original geometry to preserve, so soup_reconstruct
    // collapse is welcomed rather than rejected.
    std::vector<bool> closed_queue_preserve;
    std::vector<Mesh> manifold_closed_queue;  // bypass-Phase-8R sink (FWN+LevelSet outputs)
    std::vector<Mesh> pile;          // UNREPAIRABLE PILE
    // Parallel "synthetic" flag for pile. True when this component came from
    // FWN+LevelSet output or soup-aggregation (no user-original geometry,
    // and the rescue cascade already failed on it — so re-adding it via
    // Phase 10R.5 pile-reintegration produces marching-cubes garbage like
    // the "big cube with corner missing + internal hole" output observed
    // on soup_seed0011 et al). Pile-reintegration skips synthetic entries.
    std::vector<bool> pile_synthetic;
    bool any_manifold_ran = false;
    bool fwn_used         = false;  // see Phase 9R intersections gate

    // --- Phase 6R: aggregate ALL SOUP components into one bag, then run
    // soup_reconstruct ONCE. Per-component soup handling drops tiny shards
    // (face_count <= 3) and prevents manifold from repairing fragmented soups
    // such as 12-isolated-triangle "cube_soup" inputs. Aggregating gives the
    // CSG kernel the full geometry it needs to reconstruct a coherent solid.
    //
    // Tiny stray SOUP shards (≤ 5 faces total across all soup components) are
    // typically produced as a side-effect of duplicate-face removal on internal
    // fins / back-to-back triangles, NOT a real soup input. They are dropped
    // as noise — preserving them in the pile would spoil an otherwise-clean
    // output (regression for `multiple_fins`, `interior_face`, etc.).
    constexpr size_t kSoupNoiseThreshold = 6;
    {
        Mesh soup_bag;
        size_t soup_face_total = 0;
        std::vector<size_t> soup_comp_indices;
        for (size_t i = 0; i < components.size(); ++i) {
            if (components[i].cls != internal::ComponentClass::SOUP) continue;
            soup_comp_indices.push_back(i);
            const auto& ci = components[i];
            const uint32_t v_off = static_cast<uint32_t>(soup_bag.vertices.size());
            // Inline extract: copy referenced vertices and remap face indices.
            std::unordered_map<uint32_t, uint32_t> remap;
            remap.reserve(ci.face_indices.size() * 3);
            for (uint32_t fi : ci.face_indices) {
                const auto& f = result.mesh.faces[fi];
                std::array<uint32_t, 3> nf{};
                for (int k = 0; k < 3; ++k) {
                    auto it = remap.find(f[k]);
                    if (it == remap.end()) {
                        const uint32_t new_idx = static_cast<uint32_t>(soup_bag.vertices.size());
                        remap[f[k]] = new_idx;
                        soup_bag.vertices.push_back(result.mesh.vertices[f[k]]);
                        nf[k] = new_idx;
                    } else {
                        nf[k] = it->second;
                    }
                }
                soup_bag.faces.push_back(nf);
            }
            soup_face_total += ci.face_indices.size();
            (void)v_off; // unused; kept for potential future bag-of-bags layout
        }

        // Two reasons to drop SOUP comps without attempting reconstruction:
        //   (a) the input has OTHER content (non-SOUP components) — the SOUP
        //       parts are then almost certainly junk fragments (e.g. piercing
        //       triangles, internal fins reduced to 1 face by duplicate
        //       removal). Aggregating + failing to reconstruct would just
        //       paste them back into an otherwise-clean output.
        //   (b) total soup is below the noise threshold — small enough that
        //       it's not worth attempting CSG repair on.
        // Only when the *whole* input is SOUP and has enough faces do we run
        // soup_reconstruct. This matches the empirically-correct behavior on
        // multiple_crossing_pairs, multiple_fins, interior_face (drop) vs.
        // cube_soup, soup_seed* (aggregate).
        const bool has_non_soup = (soup_comp_indices.size() < components.size());
        if (soup_face_total > 0 && (has_non_soup || soup_face_total < kSoupNoiseThreshold)) {
            for (size_t idx : soup_comp_indices) {
                ComponentOutcome co;
                co.component_index = static_cast<int>(idx);
                co.status = ComponentStatus::dropped_zero_volume;
                co.stage  = "classify";
                co.reason = has_non_soup
                    ? "soup fragment alongside non-soup geometry"
                    : "soup shard below noise threshold";
                result.component_outcomes.push_back(co);
            }
        } else if (!soup_bag.faces.empty()) {
            if (!opts.soup_reconstruct) {
                pile.push_back(std::move(soup_bag));
                pile_synthetic.push_back(true);  // aggregated soup, no user-original
            } else {
                // Inline volume-vs-bbox test: a "real" solid has |vol| >
                // bbox_diag³ · 1e-6. Vol-0 outputs (e.g. back-to-back fin
                // pairs produced by hole-fill on isolated triangles) are
                // topologically closed but geometrically degenerate, and
                // should be treated as failures so the FWN+LevelSet
                // fallback can take a swing.
                auto has_real_volume = [](const Mesh& m) -> bool {
                    if (m.faces.empty() || m.vertices.empty()) return false;
                    auto d = internal::compute_diagnostics(m);
                    double lo[3] = { m.vertices[0][0], m.vertices[0][1], m.vertices[0][2] };
                    double hi[3] = { lo[0], lo[1], lo[2] };
                    for (const auto& v : m.vertices) {
                        for (int k = 0; k < 3; ++k) {
                            if (v[k] < lo[k]) lo[k] = v[k];
                            if (v[k] > hi[k]) hi[k] = v[k];
                        }
                    }
                    const double dx = hi[0]-lo[0], dy = hi[1]-lo[1], dz = hi[2]-lo[2];
                    const double diag_sq = dx*dx + dy*dy + dz*dz;
                    if (diag_sq <= 0.0) return false;
                    const double diag_cubed = std::pow(diag_sq, 1.5);
                    return std::abs(d.signed_volume) > diag_cubed * 1e-6;
                };

                auto sr = stages::reconstruct_soup(soup_bag);

                bool sr_kept = false;
                if (sr.was_needed && sr.success && has_real_volume(sr.mesh)) {
                    add_stage("soup_reconstruct");
                    result.notes.push_back("soup_reconstruct: aggregated " +
                        std::to_string(soup_comp_indices.size()) + " soup components (" +
                        std::to_string(soup_face_total) + " faces) -> " +
                        std::to_string(sr.faces_after) + " manifold faces");
                    closed_queue.push_back(std::move(sr.mesh));
                    closed_queue_preserve.push_back(false);  // soup-aggregation output, no user-original
                    any_manifold_ran = true;
                    sr_kept = true;
                } else if (!sr.was_needed) {
                    // Pre-processing (weld + orient) made it manifold without CSG.
                    closed_queue.push_back(std::move(sr.mesh));
                    closed_queue_preserve.push_back(false);
                    sr_kept = true;
                }

                if (!sr_kept) {
                    // Last-resort Stage 10: voxel-occupancy + manifold::LevelSet.
                    // Targets random-orientation triangle soups where the
                    // arrangement-based soup_reconstruct cannot identify an
                    // enclosed volume. Voxel-occupancy is orientation-
                    // agnostic, which is what the originals (repair_algorithms_
                    // investigation.md §1.3) prescribed for random soup —
                    // FWN gives garbage on mixed-winding triangles. Output
                    // is fully-generated MC geometry (no input triangles
                    // preserved).
                    // FWN-first, voxel-occupancy as fallback. The originals
                    // (repair_algorithms_investigation.md §1.3) prescribed
                    // voxel-occupancy for random soup; empirically we found
                    // FWN+rescue-cascade produces SMOOTHER intermediate MC
                    // output that Phase 8R's soup_reconstruct rescues into
                    // tighter final shells, giving lower Hausdorff distance
                    // on our soup_seed* fixtures. Voxel kicks in as a hard
                    // fallback for inputs where FWN fails entirely (e.g.
                    // mixed-winding inputs whose winding number is 0
                    // everywhere — for those, voxel-occupancy is orientation-
                    // agnostic and still finds enclosed pockets). Both
                    // paths share the manifold_to_welded_mesh post-treatment.
                    auto fr = stages::fwn_levelset(soup_bag);
                    if (!(fr.success && has_real_volume(fr.mesh))) {
                        fr = stages::voxel_levelset(soup_bag);
                    }
                    if (fr.success && has_real_volume(fr.mesh)) {
                        add_stage("voxel_levelset");
                        result.notes.push_back("voxel_levelset: reconstructed " +
                            std::to_string(soup_face_total) + " soup faces -> " +
                            std::to_string(fr.faces_out) + " marching-cubes faces");
                        dump("voxel_raw", fr.mesh);
                        // FWN output has float32-welded vertices but typically
                        // 50-100 genuine NM edges at marching-cubes pinch
                        // points where separately-extracted shells touch.
                        //
                        // Split into connected components and feed each to
                        // Phase 8R individually. Two reasons:
                        //   (1) Each clean (NM-free) subcomponent passes
                        //       Phase 8R directly into manifold_closed_queue.
                        //       The whole-mesh-as-one-component approach was
                        //       all-or-nothing: either the soup_reconstruct
                        //       rescue worked for everything (gain) or
                        //       failed and the whole reconstruction went to
                        //       pile (loss).
                        //   (2) Per-component rescue is more focused: a
                        //       small NM-laden chunk has a much better
                        //       chance through soup_reconstruct than a
                        //       17 000-tri reconstruction with 50 NM edges.
                        fwn_used = true;
                        auto subcomps = internal::classify_components(fr.mesh,
                            opts.soup_planarity_threshold, opts.soup_open_ratio_threshold);
                        for (const auto& sub : subcomps) {
                            if (sub.face_indices.empty()) continue;
                            closed_queue.push_back(
                                internal::extract_component(fr.mesh, sub.face_indices));
                            closed_queue_preserve.push_back(false);  // FWN synthetic output
                        }
                        any_manifold_ran = true;
                    } else {
                        result.notes.push_back("soup_reconstruct + voxel_levelset both failed on aggregated soup (" +
                            std::to_string(soup_face_total) + " faces, preserved in pile)");
                        pile.push_back(std::move(soup_bag));
                        pile_synthetic.push_back(true);  // soup-aggregation pile, no user-original
                        result.partial_failure = true;
                    }
                }
            }
        }
    }

    // --- Phase 7R: OPEN queue → nm_edge → orient → fill_holes → CLOSED or PILE ---
    for (const auto& ci : components) {
        if (ci.cls != internal::ComponentClass::OPEN) continue;
        Mesh comp = internal::extract_component(result.mesh, ci.face_indices);
        if (!opts.holes) {
            pile.push_back(std::move(comp));
            pile_synthetic.push_back(false);  // OPEN input component
            continue;
        }
        if (opts.non_manifold) { auto r = stages::fix_non_manifold_edges(comp); comp = std::move(r.mesh); }
        if (opts.orient)       { auto r = stages::orient_mesh(comp);            comp = std::move(r.mesh); }
        auto hr = stages::fill_holes(comp);
        comp = std::move(hr.mesh);
        auto diag = internal::compute_diagnostics(comp);
        if (diag.open_boundary_edges == 0) {
            if (hr.holes_filled > 0) {
                add_stage("holes");
                result.notes.push_back("holes: filled " + std::to_string(hr.holes_filled) +
                    " holes in open component (" + std::to_string(hr.faces_added) + " faces added)");
            }
            closed_queue.push_back(std::move(comp));
            closed_queue_preserve.push_back(true);  // OPEN component came from real input
        } else {
            result.notes.push_back("holes: open component could not be closed, preserved in pile");
            result.partial_failure = true;
            pile.push_back(std::move(comp));
            pile_synthetic.push_back(false);  // real OPEN input that couldn't be closed
        }
    }

    // NO_BOUNDARY components → directly to CLOSED queue.
    for (const auto& ci : components) {
        if (ci.cls != internal::ComponentClass::NO_BOUNDARY) continue;
        closed_queue.push_back(internal::extract_component(result.mesh, ci.face_indices));
        closed_queue_preserve.push_back(true);  // closed input component, real user geometry
    }

    // --- Phase 8R: CLOSED queue → orient → thin_features → nm_edge → MANIFOLD CLOSED or PILE ---
    // (manifold_closed_queue was declared earlier so FWN+LevelSet outputs can bypass Phase 8R)

    for (size_t comp_idx = 0; comp_idx < closed_queue.size(); ++comp_idx) {
        auto& comp = closed_queue[comp_idx];
        const bool preserve_original = (comp_idx < closed_queue_preserve.size())
                                       ? closed_queue_preserve[comp_idx]
                                       : true;  // default to preservation if flag missing
        if (opts.orient) { auto r = stages::orient_mesh(comp); comp = std::move(r.mesh); }

        if (opts.thin_features) {
            auto tr = stages::remove_thin_features(comp);
            if (tr.pairs_found > 0) {
                comp = std::move(tr.mesh);
                add_stage("thin_features");
                result.notes.push_back("thin_features: removed " +
                    std::to_string(tr.faces_removed) + " faces in closed component");
                auto hr2 = stages::fill_holes(comp);
                if (hr2.holes_filled > 0) { comp = std::move(hr2.mesh); add_stage("holes"); }
                if (opts.orient) { auto r = stages::orient_mesh(comp); comp = std::move(r.mesh); }
            }
        }

        if (opts.non_manifold) { auto r = stages::fix_non_manifold_edges(comp); comp = std::move(r.mesh); }

        // (Fix 2 — fill_holes between nm_edge and soup_reconstruct fallback —
        // was tried 2026-05-12. With various gates (preserve_original,
        // small-boundary thresholds) it produced small wins on simple-hole
        // fixtures but consistent regressions on intersections/spike_through
        // and similar complex-boundary cases. Net was a wash or slight loss
        // by Hausdorff median. Reverted in favour of Fix 1 alone.)
        auto diag = internal::compute_diagnostics(comp);
        if (diag.non_manifold_edges == 0 && diag.open_boundary_edges == 0) {
            manifold_closed_queue.push_back(std::move(comp));
        } else if (opts.soup_reconstruct) {
            // Catastrophic-collapse guard: soup_reconstruct is designed for
            // small fragmented input. On large NM-but-mostly-coherent meshes
            // (e.g. real_world/lobster, 124k faces) it voxelises the whole
            // thing and produces a tiny shell — losing essentially all
            // original geometry. Reject any rebuild that's < 10 % of input
            // (and < 50 faces absolute) and preserve the original in the
            // pile instead, where post-pipeline cleanup or pile-reintegration
            // gets another shot. Honors user goal: prefer preserving valid
            // input surface over volume-closing at any cost.
            // Guard fires only when (1) the component came from REAL input
            // (preserve_original flag) AND (2) input is LARGE (>= 1000 faces).
            //
            // Synthetic components (FWN+LevelSet output, soup-aggregation
            // reconstructions) skip the guard — there's no user-meaningful
            // original geometry; soup_reconstruct collapse is the working
            // rescue path for those (soup_seed* family).
            //
            // Real-input components with the guard active reject < 10 % rebuilds
            // as catastrophic collapse and preserve the original in pile
            // (lobster: 124k → 2 face collapse rejected, geometry kept).
            const size_t in_faces = comp.faces.size();
            const bool   guard_active = preserve_original && (in_faces >= 1000);
            const size_t min_keep = guard_active ? in_faces / 10 : 0;
            auto sr = stages::reconstruct_soup(comp);
            bool ok = false;
            if (sr.was_needed && sr.success && sr.mesh.faces.size() >= min_keep) {
                add_stage("soup_reconstruct");
                result.notes.push_back("soup_reconstruct: reconstructed still-NM closed component");
                Mesh r = std::move(sr.mesh);
                if (opts.orient) { auto or3 = stages::orient_mesh(r); r = std::move(or3.mesh); }
                auto d2 = internal::compute_diagnostics(r);
                if (d2.non_manifold_edges == 0 && d2.open_boundary_edges == 0) {
                    manifold_closed_queue.push_back(std::move(r));
                    any_manifold_ran = true;
                    ok = true;
                }
            } else if (guard_active && sr.was_needed && sr.success && sr.mesh.faces.size() < min_keep) {
                // Collapse detected — log so it's visible in notes.
                result.notes.push_back("soup_reconstruct: rejected catastrophic collapse (" +
                    std::to_string(in_faces) + " -> " +
                    std::to_string(sr.mesh.faces.size()) +
                    " faces, < 10 % preserved)");
            }
            if (!ok) {
                result.notes.push_back("closed component still non-manifold, preserved in pile");
                result.partial_failure = true;
                pile.push_back(std::move(comp));
                // Forward the origin flag — synthetic FWN+LevelSet output
                // that fell through Phase 8R rescue must be tagged as
                // synthetic in pile so Phase 10R.5 doesn't reintegrate it
                // back into the merged solid (it's the FWN garbage by now).
                pile_synthetic.push_back(!preserve_original);
            }
        } else {
            result.notes.push_back("closed component still non-manifold, preserved in pile");
            result.partial_failure = true;
            pile.push_back(std::move(comp));
            pile_synthetic.push_back(!preserve_original);
        }
    }

    // --- Phase 9R: Boolean Merge (MANIFOLD CLOSED queue) ---
    Mesh merged_solid;
    if (!manifold_closed_queue.empty()) {
        merged_solid = internal::merge_meshes(manifold_closed_queue);
        dump("pre_intersections", merged_solid);

        if (opts.intersections && !fwn_used) {
            // Skip when FWN+LevelSet produced the queue: its output is
            // manifold by construction; re-running CSG self-union via
            // Manifold's MeshGL (float32) introduces vertex collisions
            // that degrade the dense marching-cubes mesh into a non-
            // manifold form (empirically: 0 NM → 100s of NM on
            // soup_seed* inputs at samples_per_axis ≥ 20).
            internal::ScopedTimer t("intersections", result.stage_times_ms);
            auto pre_diag = internal::compute_diagnostics(merged_solid);
            if (pre_diag.component_count > 1 || pre_diag.non_manifold_edges > 0) {
                auto ir = stages::resolve_intersections(merged_solid);
                if (!ir.manifold_failed) {
                    auto post_diag = internal::compute_diagnostics(ir.mesh);
                    if (post_diag.non_manifold_edges <= pre_diag.non_manifold_edges) {
                        merged_solid = std::move(ir.mesh);
                        add_stage("intersections");
                        if (ir.had_intersections)
                            result.notes.push_back("intersections: manifold self-union resolved intersections (" +
                                std::to_string(ir.faces_before) + " -> " + std::to_string(ir.faces_after) + " faces)");
                        any_manifold_ran = true;
                    } else {
                        result.notes.push_back("intersections: result had more nm edges, skipped");
                    }
                } else {
                    result.notes.push_back("intersections: manifold failed on merged solid");
                }
            }
        }

        dump("post_intersections", merged_solid);

        if (opts.shells) {
            internal::ScopedTimer t("shells", result.stage_times_ms);
            auto sr = stages::analyze_shells(merged_solid);
            merged_solid = std::move(sr.mesh);
            result.component_count = sr.shells_kept;
            if (sr.shells_dropped > 0) {
                add_stage("shells");
                result.notes.push_back("shells: dropped " + std::to_string(sr.shells_dropped) +
                    " shells, kept " + std::to_string(sr.shells_kept));
            }
            for (const auto& sh : sr.shells) {
                if (sh.dropped) {
                    ComponentOutcome co;
                    co.component_index = static_cast<int>(sh.component_index);
                    co.status = sh.is_contained ? ComponentStatus::dropped_contained
                                                : ComponentStatus::dropped_zero_volume;
                    co.stage  = "shells";
                    co.reason = sh.is_contained ? "contained inside another shell" : "volume below threshold";
                    result.component_outcomes.push_back(co);
                }
            }
        }
    }

    dump("post_shells", merged_solid);

    // --- Phase 10R: Final Assembly ---
    if (!merged_solid.faces.empty()) {
        if (opts.orient)     { auto r = stages::orient_mesh(merged_solid);      merged_solid = std::move(r.mesh); }
        if (opts.degenerate) { auto r = stages::remove_degenerate(merged_solid); merged_solid = std::move(r.mesh); }
    }

    // --- Phase 10R.5: Pile reintegration via soup_reconstruct ---
    // Before appending the unrepairable pile raw, give it one more attempt
    // through the soup reconstruction pipeline (voxel oracle + T-T
    // arrangement). For real-world inputs with pile-preserved components
    // (open shells hole-fill couldn't close), the soup path can sometimes
    // close them using the same machinery that handles fragmented planar
    // soups. On success the recovered geometry is appended to merged_solid
    // as a separate component (intersections will not re-run, but having
    // a closed shell is far better than raw boundaries).
    if (!pile.empty() && opts.soup_reconstruct) {
        // Skip synthetic pile entries — FWN+LevelSet output that already
        // failed Phase 8R rescue is marching-cubes garbage and re-running
        // soup_reconstruct on it just preserves the garbage (observed on
        // soup_seed0011 et al: pile_bag had 19362 faces of "big-cube-with-
        // holes" geometry, soup_reconstruct returned ~the same, the result
        // got appended to merged_solid drowning the 426 good faces from
        // rescue cascade). Keep only real-input pile entries.
        std::vector<size_t> kept_pile_indices;
        for (size_t i = 0; i < pile.size(); ++i) {
            const bool synthetic = (i < pile_synthetic.size()) ? pile_synthetic[i] : false;
            if (!synthetic) kept_pile_indices.push_back(i);
        }
        const size_t synthetic_dropped = pile.size() - kept_pile_indices.size();
        if (synthetic_dropped > 0) {
            size_t dropped_faces = 0;
            for (size_t i = 0; i < pile.size(); ++i) {
                if (i < pile_synthetic.size() && pile_synthetic[i])
                    dropped_faces += pile[i].faces.size();
            }
            result.notes.push_back("pile-reintegration: dropped " +
                std::to_string(synthetic_dropped) + " synthetic components (" +
                std::to_string(dropped_faces) + " faces of FWN/soup garbage)");
        }
        Mesh pile_bag;
        for (size_t i : kept_pile_indices) {
            const Mesh& p = pile[i];
            const uint32_t v_off = static_cast<uint32_t>(pile_bag.vertices.size());
            for (const auto& v : p.vertices) pile_bag.vertices.push_back(v);
            for (const auto& f : p.faces)
                pile_bag.faces.push_back({f[0]+v_off, f[1]+v_off, f[2]+v_off});
        }
        // Catastrophic-collapse guard (same rule as Phase 8R). For t10k_195684
        // pile-reintegration produced ~10k faces from ~10k pile faces, so
        // typical successful runs pass this gate; the guard only fires on
        // pathological collapses like lobster. Only applies to large pile
        // bags (>= 1000 faces) — small piles are expected to consolidate
        // significantly via voxel rebuild.
        const size_t pile_faces = pile_bag.faces.size();
        const bool   pile_guard = (pile_faces >= 1000);
        const size_t pile_min   = pile_guard ? pile_faces / 10 : 0;
        auto sr = stages::reconstruct_soup(pile_bag);
        const bool collapsed = pile_guard && sr.was_needed && sr.success &&
                               sr.mesh.faces.size() < pile_min;
        if (collapsed) {
            result.notes.push_back("pile-reintegration: rejected catastrophic collapse (" +
                std::to_string(pile_faces) + " -> " +
                std::to_string(sr.mesh.faces.size()) + " faces)");
        }
        if (!collapsed &&
            ((sr.was_needed && sr.success) || (!sr.was_needed && !sr.mesh.faces.empty()))) {
            const uint32_t v_off = static_cast<uint32_t>(merged_solid.vertices.size());
            for (const auto& v : sr.mesh.vertices) merged_solid.vertices.push_back(v);
            for (const auto& f : sr.mesh.faces)
                merged_solid.faces.push_back({f[0]+v_off, f[1]+v_off, f[2]+v_off});
            add_stage("soup_reconstruct");
            result.notes.push_back("pile-reintegration: soup_reconstruct closed " +
                std::to_string(pile.size()) + " unrepairable components (" +
                std::to_string(pile_bag.faces.size()) + " -> " +
                std::to_string(sr.mesh.faces.size()) + " faces)");
            any_manifold_ran = true;
            pile.clear();  // consumed; don't raw-append below
        }
    }

    // Append remaining UNREPAIRABLE PILE (raw, no modification). Only
    // pile contents that the reintegration above didn't consume reach this
    // point. Mark as partial-failure since we couldn't close them.
    // SYNTHETIC pile entries (FWN/soup-aggregation that failed rescue) are
    // dropped even from the raw-append: they're MC garbage, not preservable
    // user geometry, and re-adding them just produces the
    // big-cube-with-holes output observed on soup_seed0011.
    if (!pile.empty()) {
        size_t synth_drop = 0;
        for (size_t i = 0; i < pile.size(); ++i) {
            const bool synthetic = (i < pile_synthetic.size()) ? pile_synthetic[i] : false;
            if (synthetic) { ++synth_drop; continue; }
            const auto& p = pile[i];
            const uint32_t v_off = static_cast<uint32_t>(merged_solid.vertices.size());
            for (const auto& v : p.vertices) merged_solid.vertices.push_back(v);
            for (const auto& f : p.faces)
                merged_solid.faces.push_back({f[0]+v_off, f[1]+v_off, f[2]+v_off});
        }
        if (synth_drop > 0) {
            result.notes.push_back("raw-pile-append: dropped " + std::to_string(synth_drop) +
                " synthetic components (MC garbage from FWN/soup that rescue couldn't fix)");
        }
        if (pile.size() > synth_drop) {
            result.partial_failure = true;
        }
    }

    result.mesh = std::move(merged_solid);
    dump("post_pile", result.mesh);

    // --- Float32 compatibility pass ---
    // Manifold-based stages (intersections, soup_reconstruct) produce vertices
    // at double-precision positions that may not be exactly representable in
    // float32.  When the STL output (float32) is read back by trimesh, nearby
    // vertices collapse and create non-manifold edges (count ≥ 3).
    // Only run when manifold actually ran.
    {
        // Trigger the float32-stability pass not only for manifold-using
        // stages (which produce double-precision vertices) but also for
        // hole-fill: a fan-fill centroid is at double precision and, while
        // we snap it to float32 at insertion time, occasional centroids can
        // still coincide with existing vertices in float32 — producing
        // valence-4 NM edges only when STL is round-tripped. The pass is
        // cheap when no float-collisions exist and no manifold rebuild is
        // needed (one snap+weld+diagnostics; loop exits on first iteration).
        bool needs_float32_pass = false;
        for (const auto& s : result.stages_applied)
            if (s == "intersections" || s == "soup_reconstruct" || s == "holes") {
                needs_float32_pass = true;
                break;
            }

        if (needs_float32_pass) {
            // Helper: snap all vertices to float32 and exact-weld by bit-exact
            // float32 position, then remove degenerate/duplicate faces.
            struct F32Key {
                uint32_t x, y, z;
                bool operator==(const F32Key& o) const {
                    return x == o.x && y == o.y && z == o.z;
                }
            };
            struct F32Hash {
                size_t operator()(const F32Key& k) const {
                    size_t h = k.x;
                    h ^= (size_t)k.y * 2654435761u + 0x9e3779b9u + (h << 6) + (h >> 2);
                    h ^= (size_t)k.z * 2246822519u + 0x9e3779b9u + (h << 6) + (h >> 2);
                    return h;
                }
            };
            auto snap_and_weld = [&](const Mesh& in) -> Mesh {
                const auto& in_v = in.vertices;
                const auto& in_f = in.faces;
                std::unordered_map<F32Key, uint32_t, F32Hash> seen;
                seen.reserve(in_v.size());
                std::vector<std::array<double, 3>> new_verts;
                new_verts.reserve(in_v.size());
                std::vector<uint32_t> remap(in_v.size());
                for (uint32_t i = 0; i < static_cast<uint32_t>(in_v.size()); ++i) {
                    float fx = static_cast<float>(in_v[i][0]);
                    float fy = static_cast<float>(in_v[i][1]);
                    float fz = static_cast<float>(in_v[i][2]);
                    uint32_t bx, by, bz;
                    std::memcpy(&bx, &fx, 4); std::memcpy(&by, &fy, 4); std::memcpy(&bz, &fz, 4);
                    if (bx == 0x80000000u) bx = 0;
                    if (by == 0x80000000u) by = 0;
                    if (bz == 0x80000000u) bz = 0;
                    F32Key key{bx, by, bz};
                    auto it = seen.find(key);
                    if (it == seen.end()) {
                        uint32_t idx = static_cast<uint32_t>(new_verts.size());
                        seen[key] = idx;
                        float cx, cy, cz;
                        std::memcpy(&cx, &bx, 4); std::memcpy(&cy, &by, 4); std::memcpy(&cz, &bz, 4);
                        new_verts.push_back({static_cast<double>(cx), static_cast<double>(cy), static_cast<double>(cz)});
                        remap[i] = idx;
                    } else {
                        remap[i] = it->second;
                    }
                }
                Mesh out;
                out.vertices = std::move(new_verts);
                out.faces.reserve(in_f.size());
                for (const auto& f : in_f) {
                    uint32_t a = remap[f[0]], b = remap[f[1]], c = remap[f[2]];
                    if (a != b && b != c && a != c)
                        out.faces.push_back({a, b, c});
                }
                return out;
            };

            // Loop: snap → degenerate → (if nm) manifold reconstruction using
            // float32 output (GetMeshGL, not GetMeshGL64).
            // Using GetMeshGL ensures output vertices are at float32-representable
            // positions, so the next snap iteration will have no collapses.
            // Bounded to 3 iterations; converges in 1-2 for typical cases.
            for (int iter = 0; iter < 3; ++iter) {
                Mesh snapped = snap_and_weld(result.mesh);
                {
                    // Sliver-aware threshold for the post-snap degenerate
                    // filter. Liepa hole-fill can produce legitimately
                    // narrow triangles whose absolute area sits below
                    // 1e-14 on large-scale meshes (e.g. t10k_71531 with
                    // mm-scale topology in a hundreds-of-mm bbox). The
                    // absolute-1e-14 default mis-classifies these as
                    // degenerate and drops them, re-opening boundaries
                    // that hole-fill had just closed. Scale-relative
                    // (~bbox² · 1e-22) catches truly-collinear zero-area
                    // tris but preserves narrow valid fills.
                    double sxsq = 0.0;
                    if (!snapped.vertices.empty()) {
                        double lo[3] = { snapped.vertices[0][0],
                                         snapped.vertices[0][1],
                                         snapped.vertices[0][2] };
                        double hi[3] = { lo[0], lo[1], lo[2] };
                        for (const auto& v : snapped.vertices) {
                            for (int k = 0; k < 3; ++k) {
                                if (v[k] < lo[k]) lo[k] = v[k];
                                if (v[k] > hi[k]) hi[k] = v[k];
                            }
                        }
                        const double dx = hi[0]-lo[0], dy = hi[1]-lo[1], dz = hi[2]-lo[2];
                        sxsq = dx*dx + dy*dy + dz*dz;
                    }
                    const double thresh = std::max(1e-30, sxsq * 1e-22);
                    auto dr = stages::remove_degenerate(snapped, thresh);
                    snapped = std::move(dr.mesh);
                }
                auto snap_diag = internal::compute_diagnostics(snapped);
                if (snap_diag.non_manifold_edges == 0 && snap_diag.open_boundary_edges == 0) {
                    result.mesh = std::move(snapped);
                    break;
                }
                // nm or boundary edges after snap: reconstruct via manifold.
                // Bypass soup_reconstruct's scale-aware pre-weld (which would
                // over-merge float32-distinct vertices). Build MeshGL directly.
                manifold::MeshGL snap_gl;
                snap_gl.numProp = 3;
                snap_gl.vertProperties.reserve(snapped.vertices.size() * 3);
                for (const auto& v : snapped.vertices) {
                    snap_gl.vertProperties.push_back(static_cast<float>(v[0]));
                    snap_gl.vertProperties.push_back(static_cast<float>(v[1]));
                    snap_gl.vertProperties.push_back(static_cast<float>(v[2]));
                }
                snap_gl.triVerts.reserve(snapped.faces.size() * 3);
                for (const auto& f : snapped.faces) {
                    snap_gl.triVerts.push_back(f[0]);
                    snap_gl.triVerts.push_back(f[1]);
                    snap_gl.triVerts.push_back(f[2]);
                }
                manifold::Manifold mf(snap_gl);
                if (mf.Status() != manifold::Manifold::Error::NoError || mf.NumTri() == 0) {
                    result.mesh = std::move(snapped);
                    break;
                }
                // Use GetMeshGL() (float32 output) so vertices are float32-exact.
                manifold::MeshGL mf_out = mf.GetMeshGL();
                if (mf_out.triVerts.empty() || mf_out.vertProperties.empty() || mf_out.numProp < 3) {
                    result.mesh = std::move(snapped);
                    break;
                }
                // Convert float32 positions to double (exact, lossless).
                Mesh rebuilt;
                const uint32_t nv = static_cast<uint32_t>(mf_out.vertProperties.size() / mf_out.numProp);
                rebuilt.vertices.reserve(nv);
                for (uint32_t vi = 0; vi < nv; ++vi) {
                    uint32_t base = vi * mf_out.numProp;
                    rebuilt.vertices.push_back({
                        static_cast<double>(mf_out.vertProperties[base + 0]),
                        static_cast<double>(mf_out.vertProperties[base + 1]),
                        static_cast<double>(mf_out.vertProperties[base + 2])
                    });
                }
                const uint32_t nt = static_cast<uint32_t>(mf_out.triVerts.size() / 3);
                rebuilt.faces.reserve(nt);
                for (uint32_t ti = 0; ti < nt; ++ti) {
                    uint32_t base = ti * 3;
                    rebuilt.faces.push_back({mf_out.triVerts[base], mf_out.triVerts[base+1], mf_out.triVerts[base+2]});
                }
                result.mesh = std::move(rebuilt);
                // Loop continues to snap the rebuilt mesh; if it's float32-clean
                // (vertices came from GetMeshGL), the next snap will be a no-op.
            }
        }
    }

    dump("post_float32", result.mesh);

    // --- Sliver collapse ---
    // Extreme-needle triangles (very low shape quality) — typically thin-strip
    // flap faces left by imperfect CSG self-intersection resolution — collapse
    // their shortest MANIFOLD edge, which removes the needle without opening a
    // boundary. A sliver whose short edge is the non-manifold edge itself is
    // skipped (collapsing it would tear the surface). Targets the residual NM
    // on inputs like t10k_1582375 (4 NM edges, each with 2 flap needles at
    // quality ≈ 0.001). Runs before the cleanup loop so any change it makes is
    // re-evaluated; gated on non_manifold so it only fires when there's
    // something to fix.
    if (opts.non_manifold && !result.mesh.faces.empty()) {
        auto pre_sl = internal::compute_diagnostics(result.mesh);
        if (pre_sl.non_manifold_edges > 0) {
            auto sl = stages::collapse_slivers(result.mesh);
            if (sl.slivers_collapsed > 0) {
                auto post_sl = internal::compute_diagnostics(sl.mesh);
                // Keep only if total defect did not rise.
                const int pre_d  = pre_sl.open_boundary_edges + pre_sl.non_manifold_edges;
                const int post_d = post_sl.open_boundary_edges + post_sl.non_manifold_edges;
                if (post_d <= pre_d) {
                    result.mesh = std::move(sl.mesh);
                    result.notes.push_back("sliver-collapse: collapsed " +
                        std::to_string(sl.slivers_collapsed) + " sliver edges (defect " +
                        std::to_string(pre_d) + " -> " + std::to_string(post_d) + ")");
                }
            }
        }
    }

    // --- Iterative topology cleanup (after float32 stabilisation) ---
    // The post-pipeline mesh, now float32-stable, may still have small
    // residual artifacts: NM edges from hole-fill fan-fan collisions, or
    // boundary edges where pile reintegration didn't fully close.
    //
    // Three independent guards keep this from regressing edge cases:
    //   1. STRICT IMPROVEMENT: an iteration must reduce bnd or nm WITHOUT
    //      increasing the other axis — protects against fill→NM cascades.
    //   2. VOLUME PRESERVATION: a non-trivial pre-iter volume cannot drop
    //      by more than half — protects against fills that close the
    //      surface into a degenerate (collapsing volume) shape, which
    //      empirically happens on inputs like `lobster` where the prior
    //      pipeline already left the mesh at a quasi-stable state.
    //   3. SNAPSHOT REVERT: if either guard fires, revert to the pre-iter
    //      mesh and exit. No partial-state output.
    // (A final tolerance-weld + degenerate/duplicate cleanup pass was tried
    // 2026-05-15 to resolve the chronic residual nm=1 on teapot /
    // plastic_vase / t10k_71531 — each is one NM edge with 4 faces that are
    // 2 real faces duplicated at coincident-but-differently-indexed
    // vertices. A bbox_diag·1e-6 weld over-merged and regressed teapot
    // nm1→nm3; a tighter weld doesn't reach the ~1e-6 duplicate spacing.
    // The right fix is a seam-aware weld that merges only across
    // component-merge seams (merge_meshes leaves duplicate boundary
    // vertices) — left for a future session. Reverted.)

    if (!result.mesh.faces.empty()) {
        for (int iter = 0; iter < 3; ++iter) {
            auto pre = internal::compute_diagnostics(result.mesh);
            if (pre.open_boundary_edges == 0 && pre.non_manifold_edges == 0) break;
            // Gate: enter the loop when there are NM edges to remove OR
            // when there's a SMALL pure-boundary residual (≤ 20 edges).
            //
            // Big pure-boundary residuals (open shells we couldn't close —
            // partial cylinders, Möbius strips, mug w/ thin features) are
            // stable end-states; re-running fill_holes on them empirically
            // introduces artifacts (nmV junctions, low-quality slivers)
            // without closing the surface.
            //
            // BUT small residuals like `bird_bath` (bnd=3), `sea_vase`
            // (bnd=3), `lobster` (bnd≈30) are inputs where the pipeline
            // got 99 % of the way there and left a tiny boundary loop
            // from fin-removal cascades. A single Liepa pass closes those
            // cleanly. The snapshot-revert guard below catches the cases
            // where it doesn't, so the worst-case is "same as before".
            //
            // Threshold raised 20 → 400 on 2026-05-15: the old limit was
            // set under the per-axis `worsened` rule, where a fill that
            // turned 32 boundary edges into 5 NM edges looked like a
            // regression and was reverted (the "mug regression at 50"
            // noted historically). The new total-defect `worsened` rule
            // (bnd+nm) correctly scores that as a 32→5 improvement, so the
            // loop can safely attempt larger boundary residuals — the
            // snapshot-revert + volume-collapse guards still catch genuine
            // regressions. 400 catches lobster (32), teapot (37), human
            // (55), mug (166); excludes only pathological huge-boundary
            // open shells (partial_cylinder etc.) which are stable
            // end-states fill_holes can't close anyway.
            const bool small_boundary_only =
                pre.non_manifold_edges == 0 &&
                pre.open_boundary_edges > 0 &&
                pre.open_boundary_edges <= 400;
            if (pre.non_manifold_edges == 0 && !small_boundary_only) break;

            // For small-boundary residuals, the boundaries are often
            // degenerate "pinch-point" triangles — vertices left over from
            // nm_vertex's split offset (~1e-6 × fan_radius). They define a
            // tiny triangle hole that Liepa rejects as degenerate. The fix
            // is to MERGE the near-coincident vertices, not FILL the loop.
            // Try a tolerance-weld first; if it eliminates the boundary,
            // we're done before fill_holes ever runs.
            if (small_boundary_only) {
                // Compute bbox for tolerance.
                double bb_lo[3] = { result.mesh.vertices[0][0],
                                    result.mesh.vertices[0][1],
                                    result.mesh.vertices[0][2] };
                double bb_hi[3] = { bb_lo[0], bb_lo[1], bb_lo[2] };
                for (const auto& v : result.mesh.vertices) {
                    for (int k = 0; k < 3; ++k) {
                        if (v[k] < bb_lo[k]) bb_lo[k] = v[k];
                        if (v[k] > bb_hi[k]) bb_hi[k] = v[k];
                    }
                }
                const double dx = bb_hi[0]-bb_lo[0];
                const double dy = bb_hi[1]-bb_lo[1];
                const double dz = bb_hi[2]-bb_lo[2];
                const double diag = std::sqrt(dx*dx + dy*dy + dz*dz);
                // Tolerance: 1e-6 of bbox diagonal. Tighter than initial
                // weld's auto-tolerance (mean_edge × 1e-4) but loose
                // enough to catch nm_vertex split residuals.
                const double wtol = diag * 1e-6;
                Mesh tw_snap = result.mesh;
                auto wr = stages::weld(result.mesh, wtol);
                result.mesh = std::move(wr.mesh);
                if (opts.degenerate) {
                    auto dr = stages::remove_degenerate(result.mesh);
                    result.mesh = std::move(dr.mesh);
                }
                auto post_w = internal::compute_diagnostics(result.mesh);
                if (post_w.open_boundary_edges < pre.open_boundary_edges &&
                    post_w.non_manifold_edges <= pre.non_manifold_edges) {
                    result.notes.push_back(
                        "cleanup-loop: tolerance-weld collapsed bnd " +
                        std::to_string(pre.open_boundary_edges) + " -> " +
                        std::to_string(post_w.open_boundary_edges));
                    add_stage("weld");
                    // continue to next iter — the loop will re-evaluate
                    continue;
                } else {
                    result.mesh = std::move(tw_snap);
                    // Fall through to fill_holes attempt.
                }
            }

            Mesh snapshot = result.mesh;
            const double pre_vol = std::abs(pre.signed_volume);

            if (opts.non_manifold && pre.non_manifold_edges > 0) {
                auto r = stages::fix_non_manifold_edges(result.mesh);
                result.mesh = std::move(r.mesh);
            }
            if (opts.orient) {
                auto r = stages::orient_mesh(result.mesh);
                result.mesh = std::move(r.mesh);
            }
            if (opts.holes) {
                auto mid = internal::compute_diagnostics(result.mesh);
                if (mid.open_boundary_edges > 0) {
                    auto r = stages::fill_holes(result.mesh);
                    result.mesh = std::move(r.mesh);
                }
            }

            auto post = internal::compute_diagnostics(result.mesh);
            // Compare TOTAL defect count (boundary + non-manifold edges),
            // not each axis independently. The old per-axis rule treated
            // ANY increase in either axis as "worsened" — so an iteration
            // that closed 13 boundary edges but introduced 1 NM edge
            // (t10k_71531) was reverted, throwing away a 13→1 net win to
            // avoid a 0→1 nm bump. Total-defect lets the loop keep that
            // trade and clean the residual NM on the next iteration
            // (where pre.nm>0 triggers fix_non_manifold_edges). The
            // iteration cap (3) and strict-improvement `!improved` break
            // still bound any fill↔nm ping-pong.
            const int pre_defect  = pre.open_boundary_edges + pre.non_manifold_edges;
            const int post_defect = post.open_boundary_edges + post.non_manifold_edges;
            const bool worsened = post_defect > pre_defect;
            const bool improved = post_defect < pre_defect;
            // Volume-preservation guard: if the pre-iter mesh had a
            // meaningful volume and the post-iter mesh lost half of it
            // (e.g. closed into a near-zero-volume degenerate shell), the
            // iteration was destructive — revert.
            const double post_vol = std::abs(post.signed_volume);
            const bool volume_collapsed = (pre_vol > 1e-9) && (post_vol < pre_vol * 0.5);

            if (worsened || volume_collapsed) {
                result.mesh = std::move(snapshot);
                break;
            }
            if (!improved) break;
        }
    }

    prof_lap("pre:T-junction");
    // --- T-junction resolution ---
    // A non-conforming edge — one face spanning an edge whose interior holds
    // collinear vertices the other side splits into sub-edges — reads as a
    // residual boundary (or non-manifold) edge that hole-fill cannot close
    // (a zero-width collinear slit is a zero-area loop). Splitting the
    // spanning face(s) through the collinear vertices pairs the sub-edges
    // 2-to-1 and closes the slit with no degenerate faces. Guarded on total
    // defect (bnd+nm) so a mis-split cannot worsen the mesh.
    if (opts.holes && !result.mesh.faces.empty()) {
        auto pre = internal::compute_diagnostics(result.mesh);
        if (pre.open_boundary_edges > 0 || pre.non_manifold_edges > 0) {
            auto tj = stages::split_tjunctions(result.mesh);
            if (tj.edges_split > 0) {
                auto post = internal::compute_diagnostics(tj.mesh);
                const int pre_d  = static_cast<int>(pre.open_boundary_edges) +
                                   static_cast<int>(pre.non_manifold_edges);
                const int post_d = static_cast<int>(post.open_boundary_edges) +
                                   static_cast<int>(post.non_manifold_edges);
                if (post_d < pre_d) {
                    result.mesh = std::move(tj.mesh);
                    add_stage("tjunction");
                    result.notes.push_back("tjunction: split " +
                        std::to_string(tj.edges_split) + " non-conforming edge(s), "
                        + std::to_string(tj.faces_split) + " sub-faces (defect " +
                        std::to_string(pre_d) + " -> " + std::to_string(post_d) + ")");
                }
            }
        }
    }

    prof_lap("pre:nm_patch_remesh");
    // --- NM-patch local remesh (Idea C) ---
    // Runs AFTER the cleanup loop: the cleanup loop's fill_holes is itself a
    // common source of residual NM edges (it closes a boundary loop but the
    // fill triangles' edges collide with nearby geometry — e.g. t10k_1582375
    // ends the cleanup loop at bnd=0 nm=4). For each surviving NM edge,
    // delete the small local patch of faces around it and re-close the
    // opened boundary loop(s) with Liepa hole-fill — replacing the structural
    // defect (doubled-surface flap / multi-sheet pinch) with a clean
    // single-layer 2-manifold disk. >99 % of input geometry is untouched.
    // Iterated up to 3× (a remesh can leave a smaller residual); kept each
    // round only if total defect strictly improves.
    if (opts.non_manifold && opts.holes && !result.mesh.faces.empty()) {
        // The residual NM on closed doubled-surface fixtures (trumpet,
        // t10k_195684, black_vase, t10k_1582375 ...) is index-space-INVISIBLE:
        // the doubled layer uses distinct vertex indices whose positions are
        // bit-identical once truncated to float32 (the STL output precision a
        // slicer sees) but differ in double precision by more than any safe
        // metric weld tolerance. A bbox·1e-7 metric weld therefore misses
        // them entirely. Snapping to float32 and welding by exact float32 bit
        // pattern fuses precisely those duplicates — no more, no less — so the
        // NM becomes visible to the index-based remesh_nm_patches below.
        if (opts.weld) {
            std::map<std::array<float, 3>, uint32_t> seen;
            std::vector<std::array<double, 3>> nv;
            nv.reserve(result.mesh.vertices.size());
            std::vector<uint32_t> remap(result.mesh.vertices.size());
            for (uint32_t i = 0; i < result.mesh.vertices.size(); ++i) {
                std::array<float, 3> key{
                    static_cast<float>(result.mesh.vertices[i][0]),
                    static_cast<float>(result.mesh.vertices[i][1]),
                    static_cast<float>(result.mesh.vertices[i][2])};
                auto it = seen.find(key);
                if (it == seen.end()) {
                    uint32_t idx = static_cast<uint32_t>(nv.size());
                    seen.emplace(key, idx);
                    nv.push_back({static_cast<double>(key[0]),
                                  static_cast<double>(key[1]),
                                  static_cast<double>(key[2])});
                    remap[i] = idx;
                } else {
                    remap[i] = it->second;
                }
            }
            Mesh welded;
            welded.vertices = std::move(nv);
            welded.faces.reserve(result.mesh.faces.size());
            for (const auto& f : result.mesh.faces) {
                uint32_t a = remap[f[0]], b = remap[f[1]], c = remap[f[2]];
                if (a != b && b != c && a != c) welded.faces.push_back({a, b, c});
            }
            result.mesh = std::move(welded);
        }
        auto defect_of = [](const Mesh& m) {
            auto d = internal::compute_diagnostics(m);
            return static_cast<int>(d.open_boundary_edges) +
                   static_cast<int>(d.non_manifold_edges);
        };
        for (int iter = 0; iter < 6; ++iter) {
            auto pre_p = internal::compute_diagnostics(result.mesh);
            if (pre_p.non_manifold_edges == 0) break;
            const int pre_d = static_cast<int>(pre_p.open_boundary_edges) +
                              static_cast<int>(pre_p.non_manifold_edges);

            // Try multiple patch ring sizes; adopt the SMALLEST that
            // strictly improves defect, fall back to larger only if smaller
            // doesn't improve. `rings=0` keeps the patch minimal (just NM-
            // incident faces) which is the right move on a tangled "knot"
            // of multiple NM edges sharing vertices, where rings=1 over-
            // deletes and Liepa creates new NM. `rings=2` gives Liepa more
            // room when the rings=1 boundary is still locked.
            // NOTE: `coplanar_expand` mode exists on remesh_nm_patches but
            // is intentionally NOT used here — it swallows whole coplanar
            // sheets and the Liepa refill then flat-caps a 3-D corner,
            // gouging volume (t10k_1582375 lost 2.7 % volume + a visible
            // dent for a cosmetic nm4→nm1). The total-defect guard does
            // not catch geometry damage.
            stages::NmPatchResult pr;
            bool improved = false;
            int post_d_best = pre_d;
            for (int rings : {1, 0, 2, 3, 4}) {
                auto pr_try = stages::remesh_nm_patches(result.mesh, rings,
                                                       /*coplanar_expand=*/false);
                if (!pr_try.applied) continue;
                int post_d_try = defect_of(pr_try.mesh);
                if (post_d_try < post_d_best) {
                    pr = std::move(pr_try);
                    post_d_best = post_d_try;
                    improved = true;
                    break;  // smallest ring that helps — use it
                }
            }
            if (!improved) break;  // no ring size improved — stop

            add_stage("nm_patch_remesh");
            result.notes.push_back("nm_patch_remesh: removed " +
                std::to_string(pr.patch_removed) + " patch faces, filled " +
                std::to_string(pr.faces_filled) + " (defect " +
                std::to_string(pre_d) + " -> " + std::to_string(post_d_best) + ")");
            result.mesh = std::move(pr.mesh);
        }
    }

    prof_lap("pre:collapse_nm");
    // --- Progressive guarded edge-collapse of residual NM flaps ---
    // Residual non-manifold geometry that patch-remesh cannot resolve (a
    // doubled near-coplanar membrane — a CSG-union seam, e.g. the trumpet
    // bell rim) is *erased* by progressively collapsing its edges. Every
    // collapse is guarded so the NM count and the open-boundary count can
    // only fall or hold — the mesh stays watertight at every step and only
    // an internal degenerate pocket is removed. Self-guarding: nm_after is
    // always ≤ nm_before, so it is kept whenever it strictly improved.
    if (opts.non_manifold && !result.mesh.faces.empty()) {
        auto pre = internal::compute_diagnostics(result.mesh);
        if (pre.non_manifold_edges > 0) {
            auto cr = stages::collapse_nm_region(result.mesh);
            if (cr.applied && cr.nm_after < cr.nm_before) {
                // Volume guard: erasing a degenerate (near-zero-volume)
                // pocket must barely move the signed volume. A large swing
                // means the collapses chewed into real geometry (e.g. welding
                // a thin wall) — reject the whole stage in that case.
                auto post = internal::compute_diagnostics(cr.mesh);
                double v0 = std::abs(pre.signed_volume);
                double dv = std::abs(post.signed_volume - pre.signed_volume);
                if (v0 < 1e-300 || dv <= 0.005 * v0) {
                    result.mesh = std::move(cr.mesh);
                    add_stage("collapse_nm");
                    result.notes.push_back("collapse_nm: " +
                        std::to_string(cr.collapses) + " collapses, nm " +
                        std::to_string(cr.nm_before) + " -> " +
                        std::to_string(cr.nm_after));
                } else {
                    result.notes.push_back("collapse_nm: rejected (volume "
                        "change " + std::to_string(dv / v0 * 100.0) + "%)");
                }
            }
        }
    }

    prof_lap("pre:nm_patch_2nd");
    // --- Late nm_patch_remesh: 2nd pass after collapse_nm ---
    // The early nm_patch_remesh runs before collapse_nm, so it sees the
    // pre-collapse state of the mesh. collapse_nm then tidies it. If
    // residual NM survives both, a 2nd nm_patch_remesh pass can sometimes
    // close it — the geometry collapse_nm produces is simpler and more
    // amenable to small-ring patching. Same ring-retry + total-defect guard.
    if (opts.holes && opts.non_manifold && !result.mesh.faces.empty()) {
        auto defect_now = [](const Mesh& m) {
            auto d = internal::compute_diagnostics(m);
            return static_cast<int>(d.open_boundary_edges) +
                   static_cast<int>(d.non_manifold_edges);
        };
        for (int iter = 0; iter < 3; ++iter) {
            auto pre_p = internal::compute_diagnostics(result.mesh);
            if (pre_p.non_manifold_edges == 0) break;
            const int pre_d = static_cast<int>(pre_p.open_boundary_edges) +
                              static_cast<int>(pre_p.non_manifold_edges);
            stages::NmPatchResult pr;
            bool improved = false;
            int post_d_best = pre_d;
            for (int rings : {1, 0, 2, 3, 4}) {
                auto pr_try = stages::remesh_nm_patches(result.mesh, rings,
                                                       /*coplanar_expand=*/false);
                if (!pr_try.applied) continue;
                int post_d_try = defect_now(pr_try.mesh);
                if (post_d_try < post_d_best) {
                    pr = std::move(pr_try); post_d_best = post_d_try;
                    improved = true; break;
                }
            }
            if (!improved) break;
            add_stage("nm_patch_remesh");
            result.notes.push_back("nm_patch_remesh (2nd): removed " +
                std::to_string(pr.patch_removed) + " patch faces, filled " +
                std::to_string(pr.faces_filled) + " (defect " +
                std::to_string(pre_d) + " -> " + std::to_string(post_d_best) + ")");
            result.mesh = std::move(pr.mesh);
        }
    }

    prof_lap("pre:nm_local_repair");
    // --- NM-local repair: proximity weld + strict back-to-back dedup ---
    // Targets residual NM patterns the wider stages leave behind:
    //   1. CSG-corner near-coincident duplicate vertices (t10k_1582375):
    //      two vertices at the same chamfer corner separated by ~0.007% of
    //      bbox, beneath any sensible global weld tolerance. NM-local
    //      proximity weld at bbox*1e-4 catches them.
    //   2. Exact-vertex-set back-to-back face pairs (bumpy_white) the wider
    //      remove_thin_features centroid+normal detector matches together
    //      with nearby non-pair faces, worsening defect and getting
    //      rejected. Strict frozenset-vertex-set match isolates the true
    //      pairs.
    // Whole-stage guard: total defect (bnd+nm) must strictly improve AND
    // signed volume must move ≤0.5% — same envelope as collapse_nm.
    auto _defect_local = [](const Mesh& m) {
        auto d = internal::compute_diagnostics(m);
        return static_cast<int>(d.open_boundary_edges) +
               static_cast<int>(d.non_manifold_edges);
    };
    if (!result.mesh.faces.empty()) {
        auto pre_d = _defect_local(result.mesh);
        if (pre_d > 0) {
            auto pre = internal::compute_diagnostics(result.mesh);
            auto nr = stages::nm_local_repair(result.mesh);
            if (nr.merges > 0 || nr.pairs_removed > 0) {
                auto post_d = _defect_local(nr.mesh);
                auto post = internal::compute_diagnostics(nr.mesh);
                const double v0 = std::abs(pre.signed_volume);
                const double dv = std::abs(post.signed_volume - pre.signed_volume);
                if (post_d < pre_d && (v0 < 1e-300 || dv <= 0.005 * v0)) {
                    add_stage("nm_local_repair");
                    result.notes.push_back("nm_local_repair: " +
                        std::to_string(nr.merges) + " merges, " +
                        std::to_string(nr.pairs_removed) +
                        " back-to-back pair(s) removed (defect " +
                        std::to_string(pre_d) + " -> " +
                        std::to_string(post_d) + ")");
                    result.mesh = std::move(nr.mesh);
                }
            }
        }
    }

    prof_lap("pre:late_thin");
    // --- Late back-to-back duplicate-face removal (centroid+normal) ---
    // Catches back-to-back pairs introduced by
    auto _defect = [](const Mesh& m) {
        auto d = internal::compute_diagnostics(m);
        return static_cast<int>(d.open_boundary_edges) +
               static_cast<int>(d.non_manifold_edges);
    };
    if (!result.mesh.faces.empty()) {
        auto pre_d = _defect(result.mesh);
        if (pre_d > 0) {
            auto tr = stages::remove_thin_features(result.mesh);
            if (tr.pairs_found > 0) {
                auto post_d = _defect(tr.mesh);
                if (post_d < pre_d) {
                    add_stage("thin_features");
                    result.notes.push_back("thin_features (late): removed " +
                        std::to_string(tr.faces_removed) +
                        " back-to-back faces (defect " +
                        std::to_string(pre_d) + " -> " +
                        std::to_string(post_d) + ")");
                    result.mesh = std::move(tr.mesh);
                }
            }
        }
    }

    prof_lap("pre:junk_drop");
    // --- Final junk-component drop ---
    // Tiny near-zero-volume connected components — back-to-back duplicate-
    // triangle "flaps" left at non-manifold seams (visually confirmed in
    // MeshMixer as separate ~4-triangle red-edged components on
    // plastic_vase / teapot) — are degenerate junk, not user geometry.
    //
    // Discriminator: a flap/degenerate-patch has near-zero volume
    // *relative to the main body* (its faces nearly cancel); a legitimate
    // small solid has a real volume. Drop components that are BOTH tiny
    // (≤ 30 faces) AND have |vol| < 1e-5 · (largest component's |vol|).
    // teapot's junk patch is 24 faces at vol 4.8e-5 vs the body's 24.8
    // (ratio 1.9e-6 — far below 1e-5); a legit detail at 1e-5 of body
    // volume is a ~2 %-of-bbox feature, which with only ≤30 faces would
    // be implausibly coarse — safe to treat the small+negligible-volume
    // combination as junk. Edge-connectivity BFS so a flap joined to the
    // body by only a single vertex is still isolated and dropped.
    if (!result.mesh.faces.empty()) {
        const uint32_t nf = static_cast<uint32_t>(result.mesh.faces.size());
        // edge → incident faces
        std::unordered_map<uint64_t, std::vector<uint32_t>> e2f;
        e2f.reserve(nf * 3);
        auto ekey = [](uint32_t a, uint32_t b) -> uint64_t {
            if (a > b) std::swap(a, b);
            return (static_cast<uint64_t>(a) << 32) | b;
        };
        for (uint32_t fi = 0; fi < nf; ++fi) {
            const auto& f = result.mesh.faces[fi];
            e2f[ekey(f[0], f[1])].push_back(fi);
            e2f[ekey(f[1], f[2])].push_back(fi);
            e2f[ekey(f[2], f[0])].push_back(fi);
        }
        std::vector<int32_t> comp(nf, -1);
        int32_t ncomp = 0;
        for (uint32_t seed = 0; seed < nf; ++seed) {
            if (comp[seed] != -1) continue;
            const int32_t c = ncomp++;
            std::vector<uint32_t> stack{seed};
            comp[seed] = c;
            while (!stack.empty()) {
                const uint32_t cur = stack.back(); stack.pop_back();
                const auto& f = result.mesh.faces[cur];
                for (int k = 0; k < 3; ++k) {
                    for (uint32_t nb : e2f[ekey(f[k], f[(k+1)%3])]) {
                        if (comp[nb] == -1) { comp[nb] = c; stack.push_back(nb); }
                    }
                }
            }
        }
        if (ncomp > 1) {
            std::vector<uint32_t> cface(ncomp, 0);
            std::vector<double>   cvol(ncomp, 0.0);
            for (uint32_t fi = 0; fi < nf; ++fi) {
                const auto& f = result.mesh.faces[fi];
                cface[comp[fi]]++;
                // signed tet volume of (origin, a, b, c) = a · (b × c) / 6
                const auto& a = result.mesh.vertices[f[0]];
                const auto& b = result.mesh.vertices[f[1]];
                const auto& c = result.mesh.vertices[f[2]];
                const double cx = b[1]*c[2] - b[2]*c[1];
                const double cy = b[2]*c[0] - b[0]*c[2];
                const double cz = b[0]*c[1] - b[1]*c[0];
                cvol[comp[fi]] += (a[0]*cx + a[1]*cy + a[2]*cz) / 6.0;
            }
            double max_vol = 0.0;
            int32_t largest_comp = 0;
            for (int32_t c = 0; c < ncomp; ++c) {
                if (std::abs(cvol[c]) > max_vol) { max_vol = std::abs(cvol[c]); }
            }
            { // largest by face count = the main body
                uint32_t mf = 0;
                for (int32_t c = 0; c < ncomp; ++c)
                    if (cface[c] > mf) { mf = cface[c]; largest_comp = c; }
            }
            const double vol_eps = max_vol * 1e-5;
            std::vector<bool> drop_comp(ncomp, false);
            uint32_t dropped_faces = 0, dropped_comps = 0;

            // Rule 1: tiny near-zero-volume flaps (teapot/plastic_vase seam junk).
            for (int32_t c = 0; c < ncomp; ++c) {
                if (cface[c] <= 30 && std::abs(cvol[c]) < vol_eps) {
                    drop_comp[c] = true;
                    dropped_faces += cface[c];
                    ++dropped_comps;
                }
            }

            // Rule 2: small components fully CONTAINED inside the main body
            // are CSG-union artifacts (visually confirmed on human: 3 small
            // internal components, each a leftover from shell unions). A
            // legitimate hollow chamber is a LARGE shell (comparable face
            // count to the body) and is excluded by the size gate; a real
            // external detail is not contained. Containment via a single
            // +x ray cast against the body's faces (odd crossings = inside).
            // Gate: ≤ 1 % of total faces AND |vol| < 1 % of body volume —
            // a real internal chamber has substantial volume, an artifact
            // ~0.
            if (ncomp > 1) {
                const uint32_t size_gate = std::max<uint32_t>(200, nf / 100);
                // ray-triangle (Möller) along +x
                auto ray_hit = [&](const std::array<double,3>& o,
                                   const std::array<double,3>& v0,
                                   const std::array<double,3>& v1,
                                   const std::array<double,3>& v2) -> bool {
                    const double e1x=v1[0]-v0[0], e1y=v1[1]-v0[1], e1z=v1[2]-v0[2];
                    const double e2x=v2[0]-v0[0], e2y=v2[1]-v0[1], e2z=v2[2]-v0[2];
                    // dir = (1,0,0); h = dir × e2
                    const double hx=0.0, hy=-e2z, hz=e2y;
                    const double a = e1x*hx + e1y*hy + e1z*hz;
                    if (a > -1e-12 && a < 1e-12) return false;
                    const double f = 1.0/a;
                    const double sx=o[0]-v0[0], sy=o[1]-v0[1], sz=o[2]-v0[2];
                    const double u = f*(sx*hx + sy*hy + sz*hz);
                    if (u < 0.0 || u > 1.0) return false;
                    // q = s × e1
                    const double qx=sy*e1z-sz*e1y, qy=sz*e1x-sx*e1z, qz=sx*e1y-sy*e1x;
                    const double vv = f*qx;  // dir·q = qx
                    if (vv < 0.0 || u+vv > 1.0) return false;
                    const double t = f*(e2x*qx + e2y*qy + e2z*qz);
                    return t > 1e-9;
                };
                for (int32_t c = 0; c < ncomp; ++c) {
                    if (c == largest_comp || drop_comp[c]) continue;
                    if (cface[c] > size_gate) continue;
                    if (std::abs(cvol[c]) >= max_vol * 0.01) continue;
                    // origin = centroid of first face of component c
                    std::array<double,3> origin{0,0,0};
                    for (uint32_t fi = 0; fi < nf; ++fi) {
                        if (comp[fi] != c) continue;
                        const auto& f = result.mesh.faces[fi];
                        for (int k = 0; k < 3; ++k)
                            for (int d = 0; d < 3; ++d)
                                origin[d] += result.mesh.vertices[f[k]][d] / 3.0;
                        break;
                    }
                    int crossings = 0;
                    for (uint32_t fi = 0; fi < nf; ++fi) {
                        if (comp[fi] != largest_comp) continue;
                        const auto& f = result.mesh.faces[fi];
                        if (ray_hit(origin, result.mesh.vertices[f[0]],
                                    result.mesh.vertices[f[1]],
                                    result.mesh.vertices[f[2]]))
                            ++crossings;
                    }
                    if (crossings % 2 == 1) {  // inside the body → artifact
                        drop_comp[c] = true;
                        dropped_faces += cface[c];
                        ++dropped_comps;
                    }
                }
            }
            if (dropped_comps > 0 && dropped_faces < nf) {
                // Drop ONLY the junk-component faces. Do NOT run
                // remove_degenerate here — it would also strip degenerate
                // faces from the MAIN body, opening boundaries that no
                // later stage fills (observed: teapot nm1 → bnd33). Unused
                // vertices left behind are harmless (compute_diagnostics
                // and write_stl both ignore unreferenced vertices).
                Mesh kept;
                kept.vertices = result.mesh.vertices;
                kept.faces.reserve(nf - dropped_faces);
                for (uint32_t fi = 0; fi < nf; ++fi)
                    if (!drop_comp[comp[fi]])
                        kept.faces.push_back(result.mesh.faces[fi]);
                result.mesh = std::move(kept);
                result.notes.push_back("final cleanup: dropped " +
                    std::to_string(dropped_comps) + " tiny zero-volume junk component(s), " +
                    std::to_string(dropped_faces) + " faces");
            }
        }
    }

    prof_lap("pre:nm_carve_refill");
    // --- NM carve + recursive pipeline re-entry ---
    // For stubborn residual NM edges that all the local stages failed to
    // resolve, carve the NM-incident faces + 1-ring halo, then RE-ENTER
    // repair() recursively on the carved mesh. The Python prototype
    // (`carve_minimal.py`) proved this: carve 41 faces on black_vase,
    // run meshseal CLI fresh on the result → nm 1 → 0 CLEAN, vol +0.03 %.
    // The full pipeline (weld → orient → nm_vertex split → … → late
    // cleanup chain) settles the carved mesh in ways a single-stage
    // fill_holes + collapse_nm chain cannot replicate.
    //
    // Gated:
    //   - opts.allow_carve_refill must be true (the recursive call sets
    //     it false, preventing infinite recursion).
    //   - recursion_depth must be < 2 (paranoia cap on top of the flag).
    //   - the in-memory mesh must have a position-visible NM edge after
    //     float32-snap-weld (the residual NM on doubled-surface fixtures
    //     is often index-space invisible until snapped to float32).
    //   - the recursive result must strictly improve total defect AND
    //     stay within 5 % volume of the pre-carve mesh.
    if (opts.allow_carve_refill && opts.recursion_depth < 2 &&
        !result.mesh.faces.empty()) {
        // Float32-snap-weld a probe copy to expose index-invisible NM.
        Mesh probe = result.mesh;
        {
            std::map<std::array<float, 3>, uint32_t> seen;
            std::vector<std::array<double, 3>> nv;
            nv.reserve(probe.vertices.size());
            std::vector<uint32_t> remap(probe.vertices.size());
            for (uint32_t i = 0; i < probe.vertices.size(); ++i) {
                std::array<float, 3> key{
                    static_cast<float>(probe.vertices[i][0]),
                    static_cast<float>(probe.vertices[i][1]),
                    static_cast<float>(probe.vertices[i][2])};
                auto it = seen.find(key);
                if (it == seen.end()) {
                    uint32_t idx = static_cast<uint32_t>(nv.size());
                    seen.emplace(key, idx); remap[i] = idx;
                    nv.push_back({static_cast<double>(key[0]),
                                  static_cast<double>(key[1]),
                                  static_cast<double>(key[2])});
                } else { remap[i] = it->second; }
            }
            Mesh welded;
            welded.vertices = std::move(nv);
            welded.faces.reserve(probe.faces.size());
            for (const auto& f : probe.faces) {
                uint32_t a = remap[f[0]], b = remap[f[1]], c = remap[f[2]];
                if (a != b && b != c && a != c) welded.faces.push_back({a,b,c});
            }
            probe = std::move(welded);
        }
        auto pre = internal::compute_diagnostics(probe);
        // Carve+refill is the right tool for SMALL residual NM (1-10 ish):
        // the carve is local, the Liepa refill knits a clean patch, the
        // recursive repair settles it. For LARGE residuals (kytka1: nm=101
        // after Phase-1 stall-break) the Liepa refill makes things WORSE —
        // it generates many new NM edges as its fill triangles cross
        // surrounding geometry. The internal collapse_nm then spends minutes
        // grinding through a much larger NM set than the input, only to be
        // rejected by the outer guard. Skip when the residual is too large
        // for carve+refill to help in any reasonable time.
        constexpr int kCarveRefillMaxNm = 20;
        if (pre.non_manifold_edges > 0 &&
            pre.non_manifold_edges <= kCarveRefillMaxNm) {
            // Carve NM-incident + 1-ring halo. We use the probe (float32-
            // welded) so that the carve targets the position-visible NM.
            auto cr = stages::nm_carve_refill(probe, /*halo_rings=*/1);
            if (cr.applied && cr.faces_carved > 0) {
                // Build a mesh consisting of the carved-and-Liepa-filled
                // result, then re-enter repair() on it.
                RepairOptions nested = opts;
                nested.recursion_depth   = opts.recursion_depth + 1;
                nested.allow_carve_refill = false;   // no further recursion
                RepairResult inner = repair(cr.mesh, nested);
                auto post = internal::compute_diagnostics(inner.mesh);
                const double v0 = std::abs(pre.signed_volume);
                const double dv = std::abs(post.signed_volume - pre.signed_volume);
                // Volume guard: reject a recursion that DESTROYED the model
                // (collapsed or ballooned). Dual scale, pass if EITHER holds:
                //   (a) dv <= 20 % of signed volume  — for normal solids;
                //   (b) dv <= 10 % of bounding-box volume — the escape hatch
                //       for thin / near-flat / heavily-defective parts,
                //       where signed volume is a tiny, unreliable scale.
                //       A long thin gnomon blade or a flat eye piece can
                //       legitimately shift a large *fraction* of its small
                //       signed volume while moving very little real space
                //       (Gnomon: 24 % of vol but only 9 % of bbox-volume;
                //       eyes3: 16 % of vol but 1.8 % of bbox-volume). The
                //       carve is local (NM faces + 1-ring, <1 % of F) so it
                //       physically cannot move much true volume regardless.
                // Genuine destruction fails both: a collapse loses volume
                // comparable to the whole box, a balloon exceeds it. The
                // result must ALSO be topologically better (nm down, bnd 0).
                double bbv = 0.0;
                if (!probe.vertices.empty()) {
                    double plo[3] = { probe.vertices[0][0], probe.vertices[0][1],
                                      probe.vertices[0][2] };
                    double phi[3] = { plo[0], plo[1], plo[2] };
                    for (const auto& pv : probe.vertices)
                        for (int k = 0; k < 3; ++k) {
                            if (pv[k] < plo[k]) plo[k] = pv[k];
                            if (pv[k] > phi[k]) phi[k] = pv[k];
                        }
                    bbv = (phi[0]-plo[0]) * (phi[1]-plo[1]) * (phi[2]-plo[2]);
                }
                const bool vol_ok = (v0 < 1e-300) ||
                                    (dv <= 0.20 * v0) ||
                                    (bbv > 0.0 && dv <= 0.10 * bbv);
                if (post.non_manifold_edges < pre.non_manifold_edges &&
                    post.open_boundary_edges == 0 && vol_ok) {
                    add_stage("nm_carve_refill");
                    result.notes.push_back("nm_carve_refill: carved " +
                        std::to_string(cr.faces_carved) +
                        " NM-incident faces, re-entered repair() recursively (nm " +
                        std::to_string(cr.nm_before) + " -> " +
                        std::to_string(post.non_manifold_edges) + ")");
                    // Adopt inner result's mesh; merge stages/notes for
                    // observability (prefix nested notes so the lineage is
                    // visible in the log).
                    result.mesh = std::move(inner.mesh);
                    for (const auto& s : inner.stages_applied)
                        if (std::find(result.stages_applied.begin(),
                                      result.stages_applied.end(), s) ==
                            result.stages_applied.end())
                            result.stages_applied.push_back(s);
                    for (const auto& n : inner.notes)
                        result.notes.push_back("[recursed] " + n);
                }
            }
        }
    }

    prof_lap("pre:destruction_fb");
    // --- Last-resort destruction fallback ---
    // The component-aware pipeline (Phase 5R-10R) assumes each connected
    // component is an independent solid-or-junk. A *fragmented* solid — one
    // object exported as many overlapping open shells (real-world multi-part
    // STL) — violates that: classification, the OPEN handler, shell filtering,
    // pile reintegration and the junk-drop each shed real geometry, and the
    // result keeps a fraction of the true volume.
    //
    // Detect that here: when the pipeline output is non-watertight, run a
    // whole-input alpha-wrap reconstruction (morphological closing — bridges
    // open seams instead of leaking through them, hugs the surface instead
    // of bloating to the hull) and adopt it ONLY when it is clean AND
    // recovers substantially more volume — i.e. the pipeline genuinely
    // destroyed the model. Geometry-preserving near-misses (small residual
    // nm, full volume intact) keep their original output because the wrap
    // volume does not exceed theirs by 1.5x. Cannot regress a CLEAN fixture:
    // the fallback never runs on watertight output.
    //
    // The 1.5x volume ratio is the destruction discriminator; a surface-
    // distance-to-input gate (alpha_wrapping_plan.md §4.4) is the planned
    // refinement — deferred because a closing, unlike voxel occupancy, does
    // not burst the bbox, so a wildly-wrong wrap is not a near-term risk.
    if (opts.soup_reconstruct && !result.mesh.faces.empty() && !mesh.faces.empty()) {
        // Check the defect state the way a slicer sees it — at float32 output
        // precision. The in-memory double mesh can be index-watertight while
        // the float32 STL output is not (distinct indices at float32-identical
        // positions). Snap+weld a copy by exact float32 bit pattern first.
        Mesh probe;
        {
            std::map<std::array<float, 3>, uint32_t> seen;
            std::vector<uint32_t> remap(result.mesh.vertices.size());
            for (uint32_t i = 0; i < result.mesh.vertices.size(); ++i) {
                std::array<float, 3> key{
                    static_cast<float>(result.mesh.vertices[i][0]),
                    static_cast<float>(result.mesh.vertices[i][1]),
                    static_cast<float>(result.mesh.vertices[i][2])};
                auto it = seen.find(key);
                if (it == seen.end()) {
                    uint32_t idx = static_cast<uint32_t>(probe.vertices.size());
                    seen.emplace(key, idx);
                    probe.vertices.push_back({static_cast<double>(key[0]),
                                              static_cast<double>(key[1]),
                                              static_cast<double>(key[2])});
                    remap[i] = idx;
                } else {
                    remap[i] = it->second;
                }
            }
            for (const auto& f : result.mesh.faces) {
                uint32_t a = remap[f[0]], b = remap[f[1]], c = remap[f[2]];
                if (a != b && b != c && a != c) probe.faces.push_back({a, b, c});
            }
        }
        auto cur = internal::compute_diagnostics(probe);
        const bool not_clean = cur.open_boundary_edges > 0 ||
                               cur.non_manifold_edges > 0;
        // Genuinely-destroyed pre-gate: only reconstruct when the pipeline
        // SHREDDED the model — left it as many disconnected components. That
        // is the signature of a fragmented-solid input the component-aware
        // pipeline could not keep whole. A geometry-preserving near-miss
        // (small residual nm, ONE component, geometry intact) is left
        // untouched — wrapping it would replace correct original geometry
        // with a coarse remesh, and alpha_wrap can bloat a near-miss enough
        // to fool a volume-only adopt test.
        const bool destroyed = not_clean && cur.component_count > 1;
        if (destroyed) {
            auto vr = stages::alpha_wrap(mesh);
            if (vr.success && !vr.mesh.faces.empty()) {
                // The marching-cubes output can carry a few non-manifold
                // edges from float32 vertex collisions in dense output — erase
                // them with the same guarded edge-collapse used elsewhere so
                // the reconstruction is a clean watertight solid.
                {
                    auto vd0 = internal::compute_diagnostics(vr.mesh);
                    if (vd0.non_manifold_edges > 0) {
                        auto vc = stages::collapse_nm_region(vr.mesh);
                        if (vc.applied && vc.nm_after < vc.nm_before)
                            vr.mesh = std::move(vc.mesh);
                    }
                }
                auto vd = internal::compute_diagnostics(vr.mesh);
                const bool vox_clean = vd.open_boundary_edges == 0 &&
                                       vd.non_manifold_edges == 0;
                const double cur_vol = std::abs(cur.signed_volume);
                const double vox_vol = std::abs(vd.signed_volume);
                if (vox_clean && vox_vol > 1.5 * cur_vol) {
                    result.mesh = std::move(vr.mesh);
                    add_stage("alpha_wrap");
                    result.notes.push_back(
                        "destruction fallback: pipeline output kept only "
                        + std::to_string(cur_vol) + " of ~"
                        + std::to_string(vox_vol) + " volume — recovered via "
                        "whole-input alpha-wrap reconstruction");
                }
            }
        }
    }

    // --- "Do no harm" guard ---
    // If the input was already a clean watertight solid but the pipeline
    // ended with a non-watertight result, the repair did net harm — restore
    // the (float32-welded) input. Real-world cause: an input that is
    // edge-manifold and closed (bnd=0, nm=0) but has non-manifold *vertices*
    // and/or several components; the nm_vertex split, weld over-merge, or the
    // manifold self-union introduce nm edges the later stages cannot fully
    // clear. Better to hand back the clean input untouched than a broken
    // mesh. Cannot regress the corpus: no fixture has a clean input that
    // meshseal currently turns non-clean and we want to keep.
    if (input_was_clean && !result.mesh.faces.empty()) {
        auto probe = float32_weld(result.mesh);
        auto pd = internal::compute_diagnostics(probe);
        if (pd.open_boundary_edges > 0 || pd.non_manifold_edges > 0) {
            result.mesh = input_welded;
            add_stage("input_restored");
            result.notes.push_back("do-no-harm: repair produced a "
                "non-watertight result from a watertight input — restored "
                "the original mesh unchanged");
        }
    }

    prof_lap("pre:final_diag");
    dump("final", result.mesh);

    // --- Final diagnostics ---
    double final_signed_vol = 0.0;
    double final_bbox_diag_sq = 0.0;
    {
        internal::ScopedTimer t("diagnostics", result.stage_times_ms);
        auto diag = internal::compute_diagnostics(result.mesh);
        result.component_count = diag.component_count;
        result.watertight  = (diag.open_boundary_edges == 0 &&
                              diag.non_manifold_edges == 0);
        result.is_volume   = result.watertight && (diag.signed_volume > 0.0);
        result.self_intersections = internal::count_self_intersections(result.mesh);
        final_signed_vol = diag.signed_volume;
        if (!result.mesh.vertices.empty()) {
            double lo[3] = { result.mesh.vertices[0][0],
                             result.mesh.vertices[0][1],
                             result.mesh.vertices[0][2] };
            double hi[3] = { lo[0], lo[1], lo[2] };
            for (const auto& v : result.mesh.vertices) {
                for (int k = 0; k < 3; ++k) {
                    if (v[k] < lo[k]) lo[k] = v[k];
                    if (v[k] > hi[k]) hi[k] = v[k];
                }
            }
            const double dx = hi[0]-lo[0], dy = hi[1]-lo[1], dz = hi[2]-lo[2];
            final_bbox_diag_sq = dx*dx + dy*dy + dz*dz;
        }
    }

    // --- Heuristic confidence score (0..1) ---
    // Multiplicative rules; each independent quality axis scales the score.
    //   * empty mesh             → 0.0 (no output)
    //   * not watertight         × 0.4 (open boundaries or NM edges)
    //   * watertight but vol≈0   × 0.5 (back-to-back fin pairs that pass
    //                                    the topology check but aren't a
    //                                    real solid)
    //   * partial_failure flag   × 0.7 (pipeline reported it couldn't
    //                                    process every input component)
    // Calibrated only loosely; treat as ordinal. A clean watertight volume
    // with no partial failure scores 1.0; a partial-but-watertight 0.7;
    // an empty mesh 0.0.
    if (result.mesh.faces.empty()) {
        result.confidence = 0.0;
    } else {
        double c = 1.0;
        if (!result.watertight) c *= 0.4;
        // Volume-relative-to-bbox check: pass if |vol| > bbox_diag³ · 1e-6.
        const double abs_vol = std::abs(final_signed_vol);
        const double bbox_diag_cubed = std::pow(final_bbox_diag_sq, 1.5);
        if (bbox_diag_cubed <= 0.0 || abs_vol < bbox_diag_cubed * 1e-6) {
            c *= 0.5;
        }
        if (result.partial_failure) c *= 0.7;
        result.confidence = c;
    }

    // --- Populate structured events from notes (backward-compat) ---
    // Each free-form note becomes a RepairEvent with severity inferred from
    // the note's content. Code is derived from the stage prefix ("weld: ",
    // "holes: ", etc.) for stable machine parsing.
    for (const auto& note : result.notes) {
        RepairEvent ev;
        const auto colon = note.find(':');
        if (colon != std::string::npos) {
            ev.stage = note.substr(0, colon);
            ev.message = (colon + 2 <= note.size())
                ? note.substr(colon + 2) : std::string{};
        } else {
            ev.stage = "";
            ev.message = note;
        }
        ev.code = ev.stage;
        ev.severity = EventSeverity::info;
        // Heuristic upgrade: notes mentioning "failed", "preserved in pile",
        // or "could not" indicate a non-fatal issue worth surfacing.
        if (ev.message.find("failed") != std::string::npos ||
            ev.message.find("preserved in pile") != std::string::npos ||
            ev.message.find("could not") != std::string::npos) {
            ev.severity = EventSeverity::warning;
        }
        result.events.push_back(std::move(ev));
    }

    return result;
}

} // namespace meshseal
