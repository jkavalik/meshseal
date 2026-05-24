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
#include "stages/coplanar_fan_drop.h"
#include "stl_io.h"
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <functional>
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

    // --- compute_diagnostics cache for result.mesh ---
    //
    // compute_diagnostics builds two std::map<,>s sized to O(F) on every
    // call (edge incidence, sorted-triple face dedup). On 787 k F that is
    // ~1-2 s per call; the pipeline issues ~12 calls on result.mesh, many
    // of them pre-stage diagnostics where the previous stage didn't touch
    // result.mesh (early return). Caching the diag and invalidating on
    // mesh mutation reclaims most of that wasted work.
    //
    // Fingerprint design: data() pointers + sizes + content samples at
    // begin/middle/end of vertices and faces. Every stage in the pipeline
    // returns a NEW Mesh (then moved into result.mesh), so the storage
    // pointer changes whenever something mutates — that alone is normally
    // enough. The content samples are a belt-and-braces guard for any
    // future in-place mutation that happens to reuse the same allocation.
    // O(1) and cheap.
    struct DiagCacheFP {
        const void* v_ptr = nullptr;
        const void* f_ptr = nullptr;
        std::size_t v_size = 0;
        std::size_t f_size = 0;
        double      v_first_x = 0.0, v_mid_x = 0.0, v_last_x = 0.0;
        uint32_t    f_first_a = 0,   f_mid_a = 0,   f_last_a = 0;
        static DiagCacheFP of(const Mesh& m) {
            DiagCacheFP fp;
            fp.v_ptr  = m.vertices.empty() ? nullptr
                       : static_cast<const void*>(m.vertices.data());
            fp.f_ptr  = m.faces.empty()    ? nullptr
                       : static_cast<const void*>(m.faces.data());
            fp.v_size = m.vertices.size();
            fp.f_size = m.faces.size();
            if (!m.vertices.empty()) {
                fp.v_first_x = m.vertices.front()[0];
                fp.v_last_x  = m.vertices.back()[0];
                fp.v_mid_x   = m.vertices[fp.v_size / 2][0];
            }
            if (!m.faces.empty()) {
                fp.f_first_a = m.faces.front()[0];
                fp.f_last_a  = m.faces.back()[0];
                fp.f_mid_a   = m.faces[fp.f_size / 2][0];
            }
            return fp;
        }
        bool operator==(const DiagCacheFP& o) const {
            return v_ptr == o.v_ptr && f_ptr == o.f_ptr
                && v_size == o.v_size && f_size == o.f_size
                && v_first_x == o.v_first_x && v_mid_x == o.v_mid_x
                && v_last_x == o.v_last_x
                && f_first_a == o.f_first_a && f_mid_a == o.f_mid_a
                && f_last_a == o.f_last_a;
        }
    };
    struct DiagCache {
        bool valid = false;
        DiagCacheFP fp;
        internal::MeshDiagnostics diag;
        std::size_t hits = 0;
        std::size_t misses = 0;
    } result_diag_cache;
    auto cached_diag = [&](const Mesh& m) -> const internal::MeshDiagnostics& {
        auto cur = DiagCacheFP::of(m);
        if (result_diag_cache.valid && result_diag_cache.fp == cur) {
            ++result_diag_cache.hits;
            return result_diag_cache.diag;
        }
        ++result_diag_cache.misses;
        result_diag_cache.diag  = internal::compute_diagnostics(m);
        result_diag_cache.fp    = cur;
        result_diag_cache.valid = true;
        return result_diag_cache.diag;
    };
    // Convenience wrapper for the common case. A typical pattern is:
    //     auto post = compute_diagnostics(stage_result.mesh);
    //     result.mesh = std::move(stage_result.mesh);
    // After the move, result.mesh holds stage_result.mesh's underlying
    // storage — same data() pointer, same fingerprint. So the NEXT pre-
    // diag on result.mesh hits the cache for free. Routing both sides
    // through cached_diag is what makes the move-then-pre-diag pattern
    // a cache hit.
    auto diag_of_result = [&]() -> const internal::MeshDiagnostics& {
        return cached_diag(result.mesh);
    };

    // --- Per-stage face-count caps for huge meshes ---
    //
    // Two manifold-CSG-based stages hang on huge inputs because
    // manifold::Manifold construction is super-linear in F (BVH build +
    // edge welding + decomposition). Cap them so huge meshes don't spend
    // hundreds of seconds in libraries we don't control.
    //
    //   reconstruct_soup (Phase 8R closed-queue + Phase 10R.5 pile):
    //     800 k F. On 1 M+ F real-input components the rebuild output
    //     usually gets rejected by the existing 10 % catastrophic-collapse
    //     guard anyway, so we just save the wasted work. Corpus largest
    //     closed-queue component is lobster ~124 k F — well below.
    //
    //   intersections (Phase 9R Boolean self-union): 1.2 M F. With the
    //     AABB-disjoint pre-check (provably a no-op for disjoint shells)
    //     this only fires for AABB-overlapping huge inputs — a real loss
    //     but no alternative within manifold's CSG.
    //
    // T-junction, cleanup_loop, nm_patch_remesh, nm_local_repair, late_thin
    // etc. are O(F) or worse but DO USEFUL WORK on huge meshes (cleanup_loop
    // on Groot's 2 M F: 47 s to close all 186 boundary edges; net output
    // CLEAN vs the cap-skipped output's bnd=186). No caps for those — pay
    // the time, get clean output.
    constexpr size_t kReconstructSoupMaxFaces = 800000;
    constexpr size_t kHugeCleanupMaxFaces     = 1200000;

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

    // --- Spatial cluster split (multi-object STL dispatcher) ---
    //
    // If the welded input decomposes into multiple AABB-disjoint clusters of
    // connected components — the signature of a multi-object STL (several
    // independent objects exported on one print plate) — recurse per-cluster.
    // Each cluster runs the full pipeline at its own bbox scale, which is
    // the meaningful win for stages whose cost or quality scales with the
    // input bbox: alpha_wrap's voxel grid is now sized per-cluster instead
    // of having to resolve every cluster inside one whole-plate grid.
    //
    // Mathematical safety: AABB-disjoint shells have zero geometric coupling.
    // The manifold-union of disjoint manifolds is concatenation, so
    // per-cluster intersections produces the same NM as the whole-mesh path.
    // No corpus regression risk: the corpus is mostly single-object meshes
    // and falls through this gate at near-zero cost.
    //
    // Gates: opts.allow_spatial_split (off in the recursive call), and a
    // minimum-cluster count of 2 (one cluster = nothing to split). The
    // recursive call gets allow_spatial_split = false to avoid re-entry.
    if (opts.allow_spatial_split && result.mesh.faces.size() > 64) {
        auto comps = internal::classify_components(
            result.mesh,
            opts.soup_planarity_threshold,
            opts.soup_open_ratio_threshold);
        if (comps.size() >= 2) {
            // Per-component AABBs (over only the vertices that component
            // references, so a junk shard far from the main body lands in
            // its own cluster correctly).
            struct AABB { double lo[3]; double hi[3]; };
            std::vector<AABB> bb(comps.size());
            for (size_t ci = 0; ci < comps.size(); ++ci) {
                AABB& b = bb[ci];
                for (int k = 0; k < 3; ++k) { b.lo[k] =  1e300; b.hi[k] = -1e300; }
                for (uint32_t fi : comps[ci].face_indices) {
                    const auto& f = result.mesh.faces[fi];
                    for (int v = 0; v < 3; ++v) {
                        const auto& p = result.mesh.vertices[f[v]];
                        for (int k = 0; k < 3; ++k) {
                            if (p[k] < b.lo[k]) b.lo[k] = p[k];
                            if (p[k] > b.hi[k]) b.hi[k] = p[k];
                        }
                    }
                }
            }
            // Global bbox for padding scale.
            double gbb_lo[3] = { bb[0].lo[0], bb[0].lo[1], bb[0].lo[2] };
            double gbb_hi[3] = { bb[0].hi[0], bb[0].hi[1], bb[0].hi[2] };
            for (size_t ci = 1; ci < bb.size(); ++ci)
                for (int k = 0; k < 3; ++k) {
                    if (bb[ci].lo[k] < gbb_lo[k]) gbb_lo[k] = bb[ci].lo[k];
                    if (bb[ci].hi[k] > gbb_hi[k]) gbb_hi[k] = bb[ci].hi[k];
                }
            const double gext = std::sqrt(
                (gbb_hi[0]-gbb_lo[0])*(gbb_hi[0]-gbb_lo[0]) +
                (gbb_hi[1]-gbb_lo[1])*(gbb_hi[1]-gbb_lo[1]) +
                (gbb_hi[2]-gbb_lo[2])*(gbb_hi[2]-gbb_lo[2]));
            // Padding: 0.1 % of global bbox diagonal. Tight enough that
            // truly separate objects stay separate, loose enough that
            // shells touching at a single edge / vertex still merge into
            // one cluster (typical CSG-style input — petals + stem of a
            // flower in kytka1).
            const double pad = gext * 1e-3;
            // Union-Find by AABB overlap-with-padding.
            std::vector<size_t> uf(comps.size());
            for (size_t i = 0; i < uf.size(); ++i) uf[i] = i;
            std::function<size_t(size_t)> find_root = [&](size_t i) -> size_t {
                while (uf[i] != i) { uf[i] = uf[uf[i]]; i = uf[i]; }
                return i;
            };
            auto aabb_overlap = [&](const AABB& a, const AABB& b2) {
                for (int k = 0; k < 3; ++k) {
                    if (a.hi[k] + pad < b2.lo[k]) return false;
                    if (b2.hi[k] + pad < a.lo[k]) return false;
                }
                return true;
            };
            for (size_t i = 0; i < bb.size(); ++i)
                for (size_t j = i + 1; j < bb.size(); ++j)
                    if (aabb_overlap(bb[i], bb[j])) {
                        size_t ri = find_root(i), rj = find_root(j);
                        if (ri != rj) uf[ri] = rj;
                    }
            // Group components by root.
            std::unordered_map<size_t, std::vector<size_t>> clusters;
            for (size_t i = 0; i < uf.size(); ++i)
                clusters[find_root(i)].push_back(i);
            // Size-aware split gate:
            //   - ≥ 3 clusters: always split (multi-object plate inputs:
            //     Marble 15 → 3 clusters, kytka1 22 → 7 — the perf and
            //     alpha_wrap-quality wins outweigh per-cluster divergence)
            //   - ≥ 2 clusters AND F > 1 M: split anyway (huge meshes,
            //     where the alternative is Phase 9R intersections taking
            //     hundreds of seconds — Groot_v1_1M_Merged at 2 M F
            //     spends 424 s in manifold construction without the split)
            //   - 2 clusters AND F ≤ 1 M: do NOT split. Main_Body_parts
            //     (271 k F, 2 AABB-disjoint components) leaves nm=5
            //     residual under per-cluster processing where the
            //     whole-mesh path gets nm=0 — small per-stage cascades.
            //     The dual-scale volume guards on collapse_nm /
            //     nm_local_repair close most of that gap but not all.
            const bool many_clusters = clusters.size() >= 3;
            const bool huge_mesh_2cluster = (clusters.size() >= 2 &&
                                             result.mesh.faces.size() > 1000000);
            if (many_clusters || huge_mesh_2cluster) {
                if (prof_env)
                    std::fprintf(stderr,
                        "[prof] spatial_split: %zu components -> %zu clusters\n",
                        comps.size(), clusters.size());
                // Extract each cluster as a sub-mesh, recurse, and combine.
                RepairOptions sub_opts = opts;
                sub_opts.allow_spatial_split = false;
                std::vector<Mesh> sub_meshes;
                bool   all_watertight = true;
                bool   all_volume     = true;
                int    total_components = 0;
                int    total_si       = 0;
                bool   any_si_unknown = false;
                bool   any_partial    = false;
                double min_confidence = 1.0;
                size_t sub_idx = 0;
                for (auto& kv : clusters) {
                    std::vector<uint32_t> face_idx;
                    for (size_t ci : kv.second)
                        face_idx.insert(face_idx.end(),
                            comps[ci].face_indices.begin(),
                            comps[ci].face_indices.end());
                    Mesh sub = internal::extract_component(result.mesh, face_idx);
                    auto sub_res = repair(sub, sub_opts);
                    if (!sub_res.mesh.faces.empty())
                        sub_meshes.push_back(std::move(sub_res.mesh));
                    if (!sub_res.watertight) all_watertight = false;
                    if (!sub_res.is_volume)  all_volume     = false;
                    total_components += sub_res.component_count;
                    if (sub_res.self_intersections < 0) any_si_unknown = true;
                    else                                 total_si += sub_res.self_intersections;
                    if (sub_res.partial_failure) any_partial = true;
                    if (sub_res.confidence < min_confidence)
                        min_confidence = sub_res.confidence;
                    for (const auto& s : sub_res.stages_applied)
                        if (std::find(result.stages_applied.begin(),
                                      result.stages_applied.end(), s) ==
                            result.stages_applied.end())
                            result.stages_applied.push_back(s);
                    for (const auto& n : sub_res.notes)
                        result.notes.push_back(
                            "[cluster " + std::to_string(sub_idx) + "] " + n);
                    ++sub_idx;
                }
                result.mesh = internal::merge_meshes(sub_meshes);
                result.notes.push_back("spatial_split: " +
                    std::to_string(comps.size()) + " input components -> " +
                    std::to_string(clusters.size()) + " AABB-disjoint clusters, "
                    "repaired independently");
                result.watertight        = all_watertight;
                result.is_volume         = all_volume;
                result.component_count   = total_components;
                result.self_intersections = any_si_unknown ? -1 : total_si;
                result.partial_failure   = any_partial;
                result.confidence        = min_confidence;
                // Skip the rest of the pipeline — sub-clusters already ran it.
                return result;
            }
        }
    }

    prof_lap("pre:orient");
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

    prof_lap("pre:nm_vertex");
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

    prof_lap("pre:nm_edge");
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
        auto diag = diag_of_result();
        result.component_count = diag.component_count;
        result.watertight = (diag.open_boundary_edges == 0 && diag.non_manifold_edges == 0);
        result.is_volume  = result.watertight && (diag.signed_volume > 0.0);
        return result;
    }

    prof_lap("pre:bridge_loops");
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

    prof_lap("pre:phase5R");
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
                auto has_real_volume = [&](const Mesh& m) -> bool {
                    if (m.faces.empty() || m.vertices.empty()) return false;
                    auto d = cached_diag(m);
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

                prof_lap("6R:pre-reconstruct_soup");
                auto sr = stages::reconstruct_soup(soup_bag);
                prof_lap("6R:post-reconstruct_soup");

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
                    // FWN cost is O(F · grid_voxels). At the default
                    // samples_per_axis=20 (8 000 query points), on a 435 k F
                    // soup bag (UMODE_Ariel: 3821 fragmented shells
                    // aggregated into one bag) that is 3.5 × 10⁹ atan2 ops
                    // ≈ 150 s.
                    //
                    // We cannot just skip FWN on these inputs — voxel_levelset's
                    // flood-fill leaks through the gaps between the 3821
                    // shells and produces vol-0 garbage. FWN's winding-number
                    // reconstruction is the only thing that actually works on
                    // fragmented multi-shell soup. Instead, drop the grid
                    // resolution adaptively so query_count × F stays bounded.
                    //
                    // Budget: ~5 × 10⁸ atan2 ops (~20 s). For F ≤ 60 000 the
                    // default samples=20 (8 000 query points) fits the
                    // budget unchanged. For larger F, cube-root of (budget/F)
                    // gives the target. Floored at 8 (8³=512 query points
                    // still produces ~1000-2000 MC faces; below that the
                    // marching-cubes surface gets too blocky to be useful).
                    int fwn_samples = 20;
                    if (soup_bag.faces.size() > 60000) {
                        const double budget = 5.0e8;
                        const double f = static_cast<double>(soup_bag.faces.size());
                        int s_int = static_cast<int>(std::cbrt(budget / f));
                        if (s_int < 8)  s_int = 8;
                        if (s_int > 20) s_int = 20;
                        fwn_samples = s_int;
                    }
                    prof_lap("6R:pre-recon");
                    auto fr = stages::fwn_levelset(soup_bag, fwn_samples);
                    bool used_fwn = (fr.success && has_real_volume(fr.mesh));
                    if (!used_fwn) {
                        fr = stages::voxel_levelset(soup_bag);
                    }
                    prof_lap("6R:post-recon");
                    if (fr.success && has_real_volume(fr.mesh)) {
                        const char* reconstructor =
                            used_fwn ? "fwn_levelset" : "voxel_levelset";
                        add_stage(reconstructor);
                        result.notes.push_back(std::string(reconstructor) +
                            ": reconstructed " +
                            std::to_string(soup_face_total) + " soup faces -> " +
                            std::to_string(fr.faces_out) + " marching-cubes faces" +
                            (fwn_samples != 20
                                ? " (samples=" + std::to_string(fwn_samples) + ")"
                                : std::string()));
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

    prof_lap("pre:phase7R");
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
        auto diag = cached_diag(comp);
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

    prof_lap("pre:phase8R");
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
        auto diag = cached_diag(comp);
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
            // Size cap on reconstruct_soup. On a still-NM closed
            // component of 1 M+ F (Groot_v1_1M_Merged has two ~1 M F
            // components), reconstruct_soup runs a voxel-oracle + the
            // legacy manifold-based CSG rebuild; the manifold construction
            // alone takes hundreds of seconds. The result for huge real-
            // input components is almost always rejected anyway (the
            // 10 %-of-input min_keep guard catches catastrophic collapse).
            // Skip the call directly: route to pile, let pile reintegration
            // or the destruction-fallback alpha_wrap take over. The cap
            // (800 k F) is well above all corpus closed-queue components
            // (largest corpus fixture is lobster ~124 k F).
            stages::SoupResult sr;
            if (in_faces > kReconstructSoupMaxFaces) {
                sr.was_needed = true;
                sr.success    = false;
            } else {
                sr = stages::reconstruct_soup(comp);
            }
            bool ok = false;
            if (sr.was_needed && sr.success && sr.mesh.faces.size() >= min_keep) {
                add_stage("soup_reconstruct");
                result.notes.push_back("soup_reconstruct: reconstructed still-NM closed component");
                Mesh r = std::move(sr.mesh);
                if (opts.orient) { auto or3 = stages::orient_mesh(r); r = std::move(or3.mesh); }
                auto d2 = cached_diag(r);
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

    prof_lap("pre:phase9R");
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
            auto pre_diag = cached_diag(merged_solid);
            // AABB-disjoint pre-check: if all merged_solid components are
            // pairwise AABB-disjoint with a small padding, the manifold
            // self-union is provably a no-op (Boolean union of disjoint
            // manifolds is just concatenation). Skip the heavyweight
            // manifold construction entirely — this is a mathematical
            // identity, not a heuristic. Cheap to check (O(C²) AABB pairs
            // where C = component count).
            bool aabb_disjoint = false;
            if (pre_diag.component_count > 1) {
                auto subs = internal::classify_components(merged_solid,
                    opts.soup_planarity_threshold,
                    opts.soup_open_ratio_threshold);
                if (subs.size() == static_cast<size_t>(pre_diag.component_count) &&
                    subs.size() >= 2 && subs.size() <= 256) {
                    struct AB { double lo[3], hi[3]; };
                    std::vector<AB> abb(subs.size());
                    for (size_t i = 0; i < subs.size(); ++i) {
                        AB& b = abb[i];
                        for (int k = 0; k < 3; ++k) { b.lo[k]= 1e300; b.hi[k]=-1e300; }
                        for (uint32_t fi : subs[i].face_indices) {
                            const auto& f = merged_solid.faces[fi];
                            for (int v = 0; v < 3; ++v) {
                                const auto& p = merged_solid.vertices[f[v]];
                                for (int k = 0; k < 3; ++k) {
                                    if (p[k] < b.lo[k]) b.lo[k] = p[k];
                                    if (p[k] > b.hi[k]) b.hi[k] = p[k];
                                }
                            }
                        }
                    }
                    // Global extent for a tiny padding.
                    double glo[3] = { abb[0].lo[0], abb[0].lo[1], abb[0].lo[2] };
                    double ghi[3] = { abb[0].hi[0], abb[0].hi[1], abb[0].hi[2] };
                    for (size_t i = 1; i < abb.size(); ++i)
                        for (int k = 0; k < 3; ++k) {
                            if (abb[i].lo[k] < glo[k]) glo[k] = abb[i].lo[k];
                            if (abb[i].hi[k] > ghi[k]) ghi[k] = abb[i].hi[k];
                        }
                    const double gext = std::sqrt(
                        (ghi[0]-glo[0])*(ghi[0]-glo[0]) +
                        (ghi[1]-glo[1])*(ghi[1]-glo[1]) +
                        (ghi[2]-glo[2])*(ghi[2]-glo[2]));
                    const double pad = gext * 1e-6;
                    aabb_disjoint = true;
                    for (size_t i = 0; i < abb.size() && aabb_disjoint; ++i)
                        for (size_t j = i+1; j < abb.size() && aabb_disjoint; ++j) {
                            bool overlap = true;
                            for (int k = 0; k < 3; ++k) {
                                if (abb[i].hi[k] + pad < abb[j].lo[k]) { overlap = false; break; }
                                if (abb[j].hi[k] + pad < abb[i].lo[k]) { overlap = false; break; }
                            }
                            if (overlap) aabb_disjoint = false;
                        }
                }
                if (aabb_disjoint) {
                    result.notes.push_back("intersections: skipped (all " +
                        std::to_string(pre_diag.component_count) +
                        " merged components are AABB-disjoint; union is a "
                        "no-op)");
                }
            }
            // Skip for huge meshes (cf. Groot_v1_1M_Merged — 424 s in
            // manifold construction at 2 M F). Trade-off: huge meshes skip
            // the SI-cleanup step, but the late-pipeline (T-junction,
            // nm_patch_remesh, collapse_nm) still cleans local NM and
            // the do-no-harm guard restores clean inputs we damage.
            const bool too_large = merged_solid.faces.size() > kHugeCleanupMaxFaces;
            if (too_large && !aabb_disjoint) {
                result.notes.push_back("intersections: skipped (mesh too "
                    "large, " + std::to_string(merged_solid.faces.size()) +
                    " faces > " + std::to_string(kHugeCleanupMaxFaces) +
                    " cap)");
            }
            if (!too_large && !aabb_disjoint &&
                (pre_diag.component_count > 1 ||
                 pre_diag.non_manifold_edges > 0)) {
                auto ir = stages::resolve_intersections(merged_solid);
                if (!ir.manifold_failed) {
                    auto post_diag = cached_diag(ir.mesh);
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

    prof_lap("pre:phase10R");
    // --- Phase 10R: Final Assembly ---
    if (!merged_solid.faces.empty()) {
        if (opts.orient)     { auto r = stages::orient_mesh(merged_solid);      merged_solid = std::move(r.mesh); }
        if (opts.degenerate) { auto r = stages::remove_degenerate(merged_solid); merged_solid = std::move(r.mesh); }
    }

    prof_lap("pre:phase10R5");
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
        // Same 800 k F cap as Phase 8R's reconstruct_soup call. A huge pile
        // is typically the leftover of Phase 8R's own size-capped skip
        // (a big NM-but-mostly-coherent closed component routed straight
        // to pile). Re-running reconstruct_soup on it here would hang in
        // exactly the same way it did in Phase 8R. Skip and let the
        // raw-pile-append + destruction-fallback paths handle it.
        stages::SoupResult sr;
        if (pile_faces > kReconstructSoupMaxFaces) {
            sr.was_needed = true;
            sr.success    = false;
            result.notes.push_back("pile-reintegration: skipped (pile too "
                "large, " + std::to_string(pile_faces) + " faces > " +
                std::to_string(kReconstructSoupMaxFaces) + " cap)");
        } else {
            sr = stages::reconstruct_soup(pile_bag);
        }
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

    prof_lap("pre:float32_compat");
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
                auto snap_diag = cached_diag(snapped);
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

    prof_lap("pre:sliver");
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
        auto pre_sl = diag_of_result();
        if (pre_sl.non_manifold_edges > 0) {
            auto sl = stages::collapse_slivers(result.mesh);
            if (sl.slivers_collapsed > 0) {
                auto post_sl = cached_diag(sl.mesh);
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

    prof_lap("pre:cleanup_loop");
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
            auto pre = diag_of_result();
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
            // Threshold raised 20 → 400 → 5000:
            //   - 20→400 (2026-05-15): old limit set under per-axis
            //     `worsened` rule that misjudged 32→5 NM as a regression.
            //     Total-defect `worsened` rule (bnd+nm) is correct.
            //   - 400→5000 (this commit): pile-reintegrated multi-shell
            //     inputs (togepi_.3mf bnd=3798, cute_ghost_parts bnd=406,
            //     test.stl bnd=408, Bee_v3.stl) had output boundary just
            //     above 400 — the cleanup_loop refused to even try.
            //     Snapshot-revert + volume-collapse guards inside the
            //     loop catch genuine regressions, and the partial_cylinder
            //     class of pathology (small bnd ~66) was never the
            //     reason for the 400 ceiling — it sits well below.
            const bool small_boundary_only =
                pre.non_manifold_edges == 0 &&
                pre.open_boundary_edges > 0 &&
                pre.open_boundary_edges <= 5000;
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
                auto post_w = diag_of_result();
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
                auto mid = diag_of_result();
                if (mid.open_boundary_edges > 0) {
                    auto r = stages::fill_holes(result.mesh);
                    result.mesh = std::move(r.mesh);
                }
            }

            auto post = diag_of_result();
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
    // T-junction is fast even on huge meshes (1.86 s on Groot's 2 M F,
    // mostly the edge-incidence map build); no face-count cap needed.
    if (opts.holes && !result.mesh.faces.empty()) {
        auto pre = diag_of_result();
        if (pre.open_boundary_edges > 0 || pre.non_manifold_edges > 0) {
            auto tj = stages::split_tjunctions(result.mesh);
            if (tj.edges_split > 0) {
                auto post = cached_diag(tj.mesh);
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
        auto defect_of = [&](const Mesh& m) {
            auto d = cached_diag(m);
            return static_cast<int>(d.open_boundary_edges) +
                   static_cast<int>(d.non_manifold_edges);
        };
        for (int iter = 0; iter < 6; ++iter) {
            auto pre_p = diag_of_result();
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
        auto pre = diag_of_result();
        if (pre.non_manifold_edges > 0) {
            auto cr = stages::collapse_nm_region(result.mesh);
            if (cr.applied && cr.nm_after < cr.nm_before) {
                // Volume guard: erasing a degenerate (near-zero-volume)
                // pocket must barely move the signed volume. A large swing
                // means the collapses chewed into real geometry (e.g. welding
                // a thin wall) — reject the whole stage in that case.
                auto post = cached_diag(cr.mesh);
                double v0 = std::abs(pre.signed_volume);
                double dv = std::abs(post.signed_volume - pre.signed_volume);
                // Dual-scale guard: accept if EITHER dv is small relative
                // to signed volume OR small relative to bbox volume.
                // Signed-volume-only fails the AABB-disjoint invariance —
                // an absolute change of (say) 700 mm³ is 0.7 % of a small
                // mesh's volume but only 0.2 % of the same mesh combined
                // with another AABB-disjoint sibling. The bbox³ secondary
                // is invariant to that splitting because each cluster's
                // bbox is unchanged whether it is processed alone or
                // alongside others. Fix discovered on Main_Body_parts
                // where spatial-split rejected a collapse the whole-mesh
                // path accepted, leaving residual nm=13 vs the whole-mesh
                // nm=0.
                double bbv0 = 0.0;
                if (!result.mesh.vertices.empty()) {
                    double lo[3] = { result.mesh.vertices[0][0],
                                     result.mesh.vertices[0][1],
                                     result.mesh.vertices[0][2] };
                    double hi[3] = { lo[0], lo[1], lo[2] };
                    for (const auto& v : result.mesh.vertices)
                        for (int k = 0; k < 3; ++k) {
                            if (v[k] < lo[k]) lo[k] = v[k];
                            if (v[k] > hi[k]) hi[k] = v[k];
                        }
                    bbv0 = (hi[0]-lo[0]) * (hi[1]-lo[1]) * (hi[2]-lo[2]);
                }
                const bool vol_ok = (v0 < 1e-300) ||
                                    (dv <= 0.005 * v0) ||
                                    (bbv0 > 0.0 && dv <= 0.005 * bbv0);
                if (vol_ok) {
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
        auto defect_now = [&](const Mesh& m) {
            auto d = cached_diag(m);
            return static_cast<int>(d.open_boundary_edges) +
                   static_cast<int>(d.non_manifold_edges);
        };
        for (int iter = 0; iter < 3; ++iter) {
            auto pre_p = diag_of_result();
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
    auto _defect_local = [&](const Mesh& m) {
        auto d = cached_diag(m);
        return static_cast<int>(d.open_boundary_edges) +
               static_cast<int>(d.non_manifold_edges);
    };
    if (!result.mesh.faces.empty()) {
        auto pre_d = _defect_local(result.mesh);
        if (pre_d > 0) {
            auto pre = diag_of_result();
            auto nr = stages::nm_local_repair(result.mesh);
            if (nr.merges > 0 || nr.pairs_removed > 0) {
                auto post_d = _defect_local(nr.mesh);
                auto post = cached_diag(nr.mesh);
                const double v0 = std::abs(pre.signed_volume);
                const double dv = std::abs(post.signed_volume - pre.signed_volume);
                // Dual-scale volume guard, same reasoning as collapse_nm:
                // signed-volume-only breaks the per-cluster invariance.
                double nbbv = 0.0;
                if (!result.mesh.vertices.empty()) {
                    double lo[3] = { result.mesh.vertices[0][0],
                                     result.mesh.vertices[0][1],
                                     result.mesh.vertices[0][2] };
                    double hi[3] = { lo[0], lo[1], lo[2] };
                    for (const auto& v : result.mesh.vertices)
                        for (int k = 0; k < 3; ++k) {
                            if (v[k] < lo[k]) lo[k] = v[k];
                            if (v[k] > hi[k]) hi[k] = v[k];
                        }
                    nbbv = (hi[0]-lo[0]) * (hi[1]-lo[1]) * (hi[2]-lo[2]);
                }
                const bool vol_ok = (v0 < 1e-300) ||
                                    (dv <= 0.005 * v0) ||
                                    (nbbv > 0.0 && dv <= 0.005 * nbbv);
                if (post_d < pre_d && vol_ok) {
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
    auto _defect = [&](const Mesh& m) {
        auto d = cached_diag(m);
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
                // Pre-compute one representative face per component (centroid
                // origin). The previous code had an inner O(F) scan inside
                // the per-candidate loop to find the first face of comp c
                // — O(C × F) total. With the cache it's O(F) once.
                std::vector<uint32_t> first_face(ncomp, UINT32_MAX);
                for (uint32_t fi = 0; fi < nf; ++fi) {
                    const int32_t c = comp[fi];
                    if (first_face[c] == UINT32_MAX) first_face[c] = fi;
                }
                // Pre-build the largest-component face index list once
                // (instead of re-scanning all F faces inside each candidate's
                // ray-cast loop). On Groot-class 2 M F meshes with many small
                // contained-candidate components this drops the containment
                // check from O(candidates × F) to O(F + candidates × F_main).
                std::vector<uint32_t> main_faces;
                main_faces.reserve(cface[largest_comp]);
                for (uint32_t fi = 0; fi < nf; ++fi)
                    if (comp[fi] == largest_comp) main_faces.push_back(fi);
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
                    // origin = centroid of pre-cached first face of comp c
                    std::array<double,3> origin{0,0,0};
                    const uint32_t ff = first_face[c];
                    if (ff != UINT32_MAX) {
                        const auto& f = result.mesh.faces[ff];
                        for (int k = 0; k < 3; ++k)
                            for (int d = 0; d < 3; ++d)
                                origin[d] += result.mesh.vertices[f[k]][d] / 3.0;
                    }
                    int crossings = 0;
                    for (uint32_t fi : main_faces) {
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

    // --- Coplanar 4-fan resolver ---
    // Drops "duplicate layer" face pairs from same-sided coplanar 4-fans
    // around NM edges (typical of multi-color 3MF exports and coplanar
    // CSG seams: 4 coplanar triangles share one chord with 2 forward +
    // 2 reverse winding — the smaller-area pair is a duplicate layer).
    // Each dropped face has 2 outer edges that connected to surrounding
    // manifold geometry, so the drop opens 2 × fans_dropped × 2 boundary
    // edges; we immediately close them with fill_holes so the net effect
    // is just "remove the NM" without leaving holes.
    // Guard: keep only if total defect (bnd + nm) strictly drops.
    if (opts.non_manifold && opts.holes && !result.mesh.faces.empty() &&
        opts.allow_carve_refill && opts.recursion_depth < 2) {
        auto cfd = stages::coplanar_fan_drop(result.mesh);
        if (cfd.applied && cfd.fans_dropped > 0) {
            // The drop leaves boundary holes that a single fill_holes
            // can't close cleanly (fills can collide with neighbour
            // surface). Re-enter repair() recursively on the dropped
            // mesh — its full late-cleanup chain (fill_holes +
            // nm_edge + orient + nm_patch_remesh + collapse_nm + …)
            // settles the result properly. Same pattern as
            // nm_carve_refill's recursive re-entry.
            auto pre = cached_diag(result.mesh);
            RepairOptions nested = opts;
            nested.recursion_depth   = opts.recursion_depth + 1;
            nested.allow_carve_refill = false;   // no further recursion
            RepairResult inner = repair(cfd.mesh, nested);
            auto post = cached_diag(inner.mesh);
            const int pre_d  = static_cast<int>(pre.open_boundary_edges) +
                               static_cast<int>(pre.non_manifold_edges);
            const int post_d = static_cast<int>(post.open_boundary_edges) +
                               static_cast<int>(post.non_manifold_edges);
            // Volume guard: dropped pairs are zero-volume (coplanar
            // duplicates), so vol should barely move. Same dual-scale
            // gate as nm_carve_refill.
            const double v0 = std::abs(pre.signed_volume);
            const double dv = std::abs(post.signed_volume - pre.signed_volume);
            double bbv = 0.0;
            if (!result.mesh.vertices.empty()) {
                double plo[3] = { result.mesh.vertices[0][0],
                                  result.mesh.vertices[0][1],
                                  result.mesh.vertices[0][2] };
                double phi[3] = { plo[0], plo[1], plo[2] };
                for (const auto& pv : result.mesh.vertices)
                    for (int k = 0; k < 3; ++k) {
                        if (pv[k] < plo[k]) plo[k] = pv[k];
                        if (pv[k] > phi[k]) phi[k] = pv[k];
                    }
                bbv = (phi[0]-plo[0]) * (phi[1]-plo[1]) * (phi[2]-plo[2]);
            }
            const bool vol_ok = (v0 < 1e-300) ||
                                (dv <= 0.05 * v0) ||
                                (bbv > 0.0 && dv <= 0.005 * bbv);
            if (post_d < pre_d && vol_ok) {
                add_stage("coplanar_fan_drop");
                result.notes.push_back("coplanar_fan_drop: dropped " +
                    std::to_string(cfd.fans_dropped) +
                    " coplanar 4-fan duplicate pair(s), recursive "
                    "cleanup (defect " + std::to_string(pre_d) + " -> " +
                    std::to_string(post_d) + ")");
                result.mesh = std::move(inner.mesh);
                for (const auto& s : inner.stages_applied)
                    if (std::find(result.stages_applied.begin(),
                                  result.stages_applied.end(), s) ==
                        result.stages_applied.end())
                        result.stages_applied.push_back(s);
                for (const auto& n : inner.notes)
                    result.notes.push_back("[cfd-recursed] " + n);
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
    // Iterate carve_refill up to 4x: each successful pass may unlock
    // further reduction. Observed: front_7color.3mf needs 2 passes
    // (nm 6→1→0), Benchy_Hull needs 3 passes (nm 3→...→0). Cap at 4 —
    // each iter runs a full recursive repair() on the carved mesh, but
    // the loop's strict-improvement check terminates early, so cases
    // that don't benefit see at most 1 wasted iter.
    int prev_iter_nm = INT_MAX;
    for (int carve_iter = 0;
         carve_iter < 4 && opts.allow_carve_refill &&
         opts.recursion_depth < 2 && !result.mesh.faces.empty();
         ++carve_iter) {
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
        auto pre = cached_diag(probe);
        // Stall break: if previous iter didn't reduce nm, stop. (Also
        // covers the nm == 0 case implicitly — INT_MAX > 0 on iter 0,
        // then after a successful pass that brought nm to 0 the loop
        // body's gate `pre.non_manifold_edges > 0` won't fire next iter.)
        if (carve_iter > 0 && (int)pre.non_manifold_edges >= prev_iter_nm)
            break;
        prev_iter_nm = (int)pre.non_manifold_edges;
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
                auto post = cached_diag(inner.mesh);
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
        auto cur = cached_diag(probe);
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
        // Early-skip when adopt criteria can't possibly fire — saves the
        // 100+ s alpha_wrap call on huge near-clean inputs (Groot 2 M F
        // with cur.nm=0 bnd>0 c>1 was eating 173 s here for a wrap that
        // wouldn't be adopted). Predictions:
        //   adopt_vol_recovery needs vox_vol > 1.5 × cur_vol. alpha_wrap
        //     preserves the input's volume, so when cur_vol ≈ input_vol
        //     (within 1.5×) vox_vol can't exceed 1.5 × cur_vol. Skip.
        //   adopt_nm_cleanup needs cur.nm > 10. Skip when cur.nm ≤ 10.
        // The early skip combines both: skip iff both predictions hold.
        double cur_vol_est = std::abs(cur.signed_volume);
        bool skip_wrap = false;
        if (destroyed && cur.non_manifold_edges <= 10) {
            // Compute input vol cheaply (we already have input_welded).
            double in_vol = 0.0;
            for (const auto& f : mesh.faces) {
                const auto& a = mesh.vertices[f[0]];
                const auto& b = mesh.vertices[f[1]];
                const auto& c = mesh.vertices[f[2]];
                const double cx = b[1]*c[2] - b[2]*c[1];
                const double cy = b[2]*c[0] - b[0]*c[2];
                const double cz = b[0]*c[1] - b[1]*c[0];
                in_vol += (a[0]*cx + a[1]*cy + a[2]*cz) / 6.0;
            }
            in_vol = std::abs(in_vol);
            if (in_vol < 1.5 * cur_vol_est) {
                skip_wrap = true;
                result.notes.push_back("destruction fallback: skipped (cur nm=" +
                    std::to_string(cur.non_manifold_edges) +
                    " <= 10 and input volume ratio " +
                    std::to_string(in_vol / std::max(cur_vol_est, 1.0)) +
                    " < 1.5; adopt criteria can't fire)");
            }
        }
        if (destroyed && !skip_wrap) {
            auto vr = stages::alpha_wrap(mesh);
            if (vr.success && !vr.mesh.faces.empty()) {
                // The marching-cubes output can carry a few non-manifold
                // edges from float32 vertex collisions in dense output — erase
                // them. collapse_nm_region handles most cases; for the
                // multi-shell-overlap inputs (kytka1: 22 components welded by
                // intersections into a single non-manifold surface) the
                // residual NM is structural and collapse_nm cannot clear it.
                // Fall through to nm_patch_remesh (BFS-delete + Liepa refill)
                // which knits a clean 2-manifold disk over each NM region.
                {
                    auto vd0 = cached_diag(vr.mesh);
                    if (vd0.non_manifold_edges > 0) {
                        auto vc = stages::collapse_nm_region(vr.mesh);
                        if (vc.applied && vc.nm_after < vc.nm_before)
                            vr.mesh = std::move(vc.mesh);
                    }
                }
                // Iterative cleanup on alpha_wrap output:
                //   nm_patch_remesh (rings=1, smallest BFS patch — clears
                //   most MC pinch-point NM cleanly; larger rings tend to
                //   CREATE new NM where Liepa fills cross neighbouring
                //   surface; observed 8 → 184 on kytka1 with rings=2)
                //   then collapse_nm (mops up residual NM after fills),
                //   repeat while total nm strictly drops.
                {
                    int prev_nm = cached_diag(vr.mesh).non_manifold_edges;
                    for (int it = 0; it < 4 && prev_nm > 0; ++it) {
                        // Ring-retry: smallest ring that strictly improves
                        // nm_after wins. Same ordering as the early
                        // nm_patch_remesh loop — empirically right (1 first,
                        // 0 next to escape tangled-NM stalls, 2-4 last resort).
                        for (int rings : {1, 0, 2, 3, 4}) {
                            auto pr = stages::remesh_nm_patches(vr.mesh, rings);
                            if (pr.applied && pr.nm_after < pr.nm_before) {
                                vr.mesh = std::move(pr.mesh);
                                break;
                            }
                        }
                        auto vc = stages::collapse_nm_region(vr.mesh);
                        if (vc.applied && vc.nm_after < vc.nm_before)
                            vr.mesh = std::move(vc.mesh);
                        const int after = cached_diag(vr.mesh).non_manifold_edges;
                        if (after >= prev_nm) break;
                        prev_nm = after;
                    }
                }
                auto vd = cached_diag(vr.mesh);
                const bool vox_clean = vd.open_boundary_edges == 0 &&
                                       vd.non_manifold_edges == 0;
                const double cur_vol = std::abs(cur.signed_volume);
                const double vox_vol = std::abs(vd.signed_volume);
                // Two adopt paths:
                //   (1) Volume-recovery: pipeline output kept only a fraction
                //       of the true volume (vox_vol > 1.5x cur_vol). Original
                //       destruction-fallback case (LDraw soup, fragmented
                //       multi-shell exports).
                //   (2) NM-cleanup: pipeline output is heavily non-manifold
                //       (cur.nm > 10), the wrap reconstruction is dramatically
                //       cleaner (≤ 1/10th the NM AND ≤ 5 absolute), AND it
                //       preserves the volume (0.7..1.5 × current). This is
                //       kytka1's pattern: 22 components welded by intersections
                //       into nm=101, wrapped to nm=4 — a 25× improvement,
                //       both still non-manifold but the wrapped output is
                //       slicer-recoverable where nm=101 is not. The 1/10
                //       ratio + absolute-5 cap mean we only adopt a wrap
                //       that's a substantial improvement; small residuals
                //       (cur.nm < 10) keep their original-geometry output.
                const bool adopt_vol_recovery = vox_clean && vox_vol > 1.5 * cur_vol;
                const bool adopt_nm_cleanup   =
                    cur.non_manifold_edges > 10 &&
                    vd.open_boundary_edges == 0 &&
                    vd.non_manifold_edges <= cur.non_manifold_edges / 10 &&
                    vd.non_manifold_edges <= 5 &&
                    vox_vol > 0.7 * cur_vol &&
                    vox_vol < 1.5 * cur_vol;
                if (adopt_vol_recovery || adopt_nm_cleanup) {
                    result.mesh = std::move(vr.mesh);
                    add_stage("alpha_wrap");
                    if (adopt_vol_recovery) {
                        result.notes.push_back(
                            "destruction fallback: pipeline output kept only "
                            + std::to_string(cur_vol) + " of ~"
                            + std::to_string(vox_vol) + " volume — recovered "
                            "via whole-input alpha-wrap reconstruction");
                    } else {
                        result.notes.push_back(
                            "destruction fallback: pipeline output had nm=" +
                            std::to_string(cur.non_manifold_edges) +
                            " — replaced with alpha-wrap reconstruction "
                            "(nm=" + std::to_string(vd.non_manifold_edges) +
                            ", volume " + std::to_string(cur_vol) + " -> " +
                            std::to_string(vox_vol) + ")");
                    }
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
        auto pd = cached_diag(probe);
        if (pd.open_boundary_edges > 0 || pd.non_manifold_edges > 0) {
            result.mesh = input_welded;
            add_stage("input_restored");
            result.notes.push_back("do-no-harm: repair produced a "
                "non-watertight result from a watertight input — restored "
                "the original mesh unchanged");
        }
    }

    prof_lap("pre:final_diag");
    if (prof_env) {
        std::fprintf(stderr, "[prof] diag_cache hits=%zu misses=%zu\n",
                     result_diag_cache.hits, result_diag_cache.misses);
    }
    dump("final", result.mesh);

    // --- Final diagnostics ---
    double final_signed_vol = 0.0;
    double final_bbox_diag_sq = 0.0;
    {
        internal::ScopedTimer t("diagnostics", result.stage_times_ms);
        auto diag = diag_of_result();
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
