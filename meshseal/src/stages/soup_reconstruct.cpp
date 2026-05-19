#include "soup.h"
#include "weld.h"
#include "orient.h"
#include "holes.h"
#include "nm_edge.h"
#include "../internal/arrangement.h"
#include "../internal/diagnostics.h"
#include "../internal/vec3.h"
#include "../internal/voxel_grid.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <manifold/manifold.h>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>
#define SOUPT(...) do { if (std::getenv("MESHSEAL_TRACE")) { std::fprintf(stderr, "[soup] "); std::fprintf(stderr, __VA_ARGS__); std::fprintf(stderr, "\n"); } } while(0)

namespace meshseal::stages {

SoupResult reconstruct_soup(const Mesh& mesh) {
    SoupResult result;
    result.faces_before = static_cast<uint32_t>(mesh.faces.size());
    result.faces_after  = result.faces_before;
    result.success      = false;
    result.was_needed   = false;

    // Step 1: check if already manifold
    auto diag = internal::compute_diagnostics(mesh);
    bool already_manifold = (diag.open_boundary_edges == 0 && diag.non_manifold_edges == 0);
    if (already_manifold) {
        result.mesh = mesh;
        return result;
    }

    result.was_needed = true;

    // ─────────────────────────────────────────────────────────────────────
    // Step 2 (NEW): voxel-oracle pre-filter.
    //
    // Build a 3-D occupancy grid, rasterize every input triangle into it,
    // 6-connected flood-fill the exterior. Each input triangle is then
    // classified geometrically: sample voxel labels at small offsets above
    // and below its centroid along its face normal. If one side is Outside
    // (reached by flood-fill) and the other is Inside (an interior pocket
    // enclosed by the input geometry), the triangle is on the boundary of
    // the inferred volume and we keep it. Otherwise it is junk — interior
    // dust or exterior stray — and we discard it.
    //
    // The grid is purely an oracle. Output triangles still come from the
    // original input; their positions are unchanged. Aggressive welding of
    // the FULL soup is avoided when filtering is enough to reach a closed
    // manifold via existing weld + orient + hole-fill stages.
    //
    // If voxel-filter + standard repair produces a watertight result, we
    // return success. Otherwise fall through to the legacy aggressive-weld
    // + manifold::Manifold path below.
    // ─────────────────────────────────────────────────────────────────────
    {
        // bbox of the input
        std::array<double, 3> g_lo = { std::numeric_limits<double>::max(),
                                       std::numeric_limits<double>::max(),
                                       std::numeric_limits<double>::max() };
        std::array<double, 3> g_hi = { std::numeric_limits<double>::lowest(),
                                       std::numeric_limits<double>::lowest(),
                                       std::numeric_limits<double>::lowest() };
        for (const auto& v : mesh.vertices) {
            for (int k = 0; k < 3; ++k) {
                g_lo[k] = std::min(g_lo[k], v[k]);
                g_hi[k] = std::max(g_hi[k], v[k]);
            }
        }
        const double g_dx = g_hi[0] - g_lo[0];
        const double g_dy = g_hi[1] - g_lo[1];
        const double g_dz = g_hi[2] - g_lo[2];
        const double g_diag = std::sqrt(g_dx * g_dx + g_dy * g_dy + g_dz * g_dz);

        // Skip voxel filter for tiny / degenerate bboxes — nothing to
        // classify and the offsets would be unstable.
        if (g_diag > 1e-9) {
            internal::VoxelGrid grid(g_lo, g_hi, /*target=*/64);

            for (const auto& f : mesh.faces) {
                grid.rasterize(mesh.vertices[f[0]],
                               mesh.vertices[f[1]],
                               mesh.vertices[f[2]]);
            }
            grid.flood_fill_exterior();

            // Sample 2 cells off the surface so we clear the conservatively-
            // marked Surface band on either side.
            const double off = grid.cell_size() * 2.0;

            Mesh filtered;
            filtered.vertices = mesh.vertices;
            filtered.faces.reserve(mesh.faces.size());

            for (const auto& f : mesh.faces) {
                const auto& a = mesh.vertices[f[0]];
                const auto& b = mesh.vertices[f[1]];
                const auto& c = mesh.vertices[f[2]];
                std::array<double, 3> cen = {
                    (a[0] + b[0] + c[0]) / 3.0,
                    (a[1] + b[1] + c[1]) / 3.0,
                    (a[2] + b[2] + c[2]) / 3.0
                };
                // Face normal
                const double e1[3] = { b[0]-a[0], b[1]-a[1], b[2]-a[2] };
                const double e2[3] = { c[0]-a[0], c[1]-a[1], c[2]-a[2] };
                double n[3] = {
                    e1[1]*e2[2] - e1[2]*e2[1],
                    e1[2]*e2[0] - e1[0]*e2[2],
                    e1[0]*e2[1] - e1[1]*e2[0]
                };
                const double nlen = std::sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);
                if (nlen < 1e-30) {
                    // Degenerate triangle — drop.
                    continue;
                }
                n[0] /= nlen; n[1] /= nlen; n[2] /= nlen;

                std::array<double, 3> p_pos = {
                    cen[0] + n[0] * off, cen[1] + n[1] * off, cen[2] + n[2] * off
                };
                std::array<double, 3> p_neg = {
                    cen[0] - n[0] * off, cen[1] - n[1] * off, cen[2] - n[2] * off
                };
                const auto lab_pos = grid.sample(p_pos);
                const auto lab_neg = grid.sample(p_neg);

                // Boundary candidate: one side Outside, the other Inside.
                // Surface-on-both-sides is ambiguous (densely packed
                // geometry) — keep conservatively so hole-fill has more
                // material to work with.
                using L = internal::VoxelGrid::Label;
                bool keep = false;
                if ((lab_pos == L::Outside && lab_neg == L::Inside) ||
                    (lab_pos == L::Inside  && lab_neg == L::Outside)) {
                    keep = true;
                } else if (lab_pos == L::Surface && lab_neg == L::Surface) {
                    // Both sides hit the Surface band — sample farther out
                    // along the normal to escape the band.
                    const double off2 = grid.cell_size() * 4.0;
                    std::array<double, 3> q_pos = {
                        cen[0] + n[0] * off2, cen[1] + n[1] * off2, cen[2] + n[2] * off2
                    };
                    std::array<double, 3> q_neg = {
                        cen[0] - n[0] * off2, cen[1] - n[1] * off2, cen[2] - n[2] * off2
                    };
                    const auto q_lp = grid.sample(q_pos);
                    const auto q_ln = grid.sample(q_neg);
                    if (q_lp != q_ln && q_lp != L::Surface && q_ln != L::Surface) {
                        keep = true;
                    }
                }

                if (keep) filtered.faces.push_back(f);
            }
            SOUPT("voxel filter: %zu -> %zu faces (cell=%g)",
                  mesh.faces.size(), filtered.faces.size(), grid.cell_size());

            // ─────────────────────────────────────────────────────────────
            // Step 2 (NEW): triangle-triangle intersection on voxel-pruned
            // candidate pairs.
            //
            // Build a cell→triangle hash; only test pairs that share at
            // least one voxel cell. For each pair that actually intersects,
            // record the intersection segment as a cut on BOTH triangles.
            // Then split each triangle by its accumulated cuts. Finally,
            // re-classify each FRAGMENT through the voxel oracle (the
            // parent's classification doesn't transfer — `cube_soup_
            // intersecting`'s face-pairs are boundary-candidates as
            // wholes, but their "extending past the cube edge" fragments
            // are entirely outside the inferred volume and must be
            // dropped here).
            // ─────────────────────────────────────────────────────────────
            const double cell = grid.cell_size();

            // Per-triangle accumulated cuts and intersection-pair tracking.
            std::vector<std::vector<std::pair<internal::Vec3d, internal::Vec3d>>>
                tri_cuts(filtered.faces.size());

            if (!filtered.faces.empty()) {
                // Hash filtered triangles into voxel cells (using each
                // triangle's bbox; cheap, conservative).
                auto cell_key = [](int i, int j, int k) -> std::uint64_t {
                    return (static_cast<std::uint64_t>(i & 0xFFFFF)) |
                           (static_cast<std::uint64_t>(j & 0xFFFFF) << 20) |
                           (static_cast<std::uint64_t>(k & 0xFFFFF) << 40);
                };
                std::unordered_map<std::uint64_t, std::vector<std::uint32_t>> cell_tris;
                cell_tris.reserve(filtered.faces.size() * 4);
                for (std::uint32_t fi = 0; fi < filtered.faces.size(); ++fi) {
                    const auto& f = filtered.faces[fi];
                    const auto& va = filtered.vertices[f[0]];
                    const auto& vb = filtered.vertices[f[1]];
                    const auto& vc = filtered.vertices[f[2]];
                    double tlo[3], thi[3];
                    for (int k = 0; k < 3; ++k) {
                        tlo[k] = std::min({va[k], vb[k], vc[k]});
                        thi[k] = std::max({va[k], vb[k], vc[k]});
                    }
                    int ilo[3], ihi[3];
                    for (int k = 0; k < 3; ++k) {
                        ilo[k] = static_cast<int>(std::floor((tlo[k] - g_lo[k]) / cell));
                        ihi[k] = static_cast<int>(std::floor((thi[k] - g_lo[k]) / cell));
                    }
                    for (int i = ilo[0]; i <= ihi[0]; ++i)
                    for (int j = ilo[1]; j <= ihi[1]; ++j)
                    for (int k = ilo[2]; k <= ihi[2]; ++k) {
                        cell_tris[cell_key(i, j, k)].push_back(fi);
                    }
                }

                // Walk candidate pairs (those sharing a cell) and compute
                // T-T intersection. Each pair seen at most once.
                std::unordered_set<std::uint64_t> seen_pairs;
                seen_pairs.reserve(filtered.faces.size() * 4);
                auto pair_key = [](std::uint32_t a, std::uint32_t b) {
                    if (a > b) std::swap(a, b);
                    return (static_cast<std::uint64_t>(a) << 32) |
                           static_cast<std::uint64_t>(b);
                };
                for (const auto& bucket : cell_tris) {
                    const auto& v = bucket.second;
                    if (v.size() < 2) continue;
                    for (size_t i = 0; i < v.size(); ++i)
                    for (size_t j = i + 1; j < v.size(); ++j) {
                        const std::uint32_t fa = v[i];
                        const std::uint32_t fb = v[j];
                        if (!seen_pairs.insert(pair_key(fa, fb)).second) continue;

                        const auto& FA = filtered.faces[fa];
                        const auto& FB = filtered.faces[fb];
                        auto seg = internal::tri_tri_intersect(
                            filtered.vertices[FA[0]], filtered.vertices[FA[1]], filtered.vertices[FA[2]],
                            filtered.vertices[FB[0]], filtered.vertices[FB[1]], filtered.vertices[FB[2]]);
                        if (seg.size() == 2) {
                            tri_cuts[fa].push_back({seg[0], seg[1]});
                            tri_cuts[fb].push_back({seg[0], seg[1]});
                        }
                    }
                }
                size_t total_cuts = 0;
                for (auto& v : tri_cuts) total_cuts += v.size();
                SOUPT("T-T intersect: %zu cut-records added across %zu filtered tris",
                      total_cuts, filtered.faces.size());

                // Cut consolidation. Möller returns the intersection clipped
                // to BOTH triangles of a pair, so a chord that physically
                // spans triangle T edge-to-edge is split across multiple
                // (T, other_i) calls — each producing a partial segment with
                // one endpoint on T's edge and one in T's interior. The
                // splitter can only handle edge-to-edge cuts. Merge cuts
                // that lie on the same line within a triangle's plane:
                // they form one logical chord and consolidate cleanly.
                const double line_tol = cell * 0.02;
                for (auto& cuts : tri_cuts) {
                    if (cuts.size() <= 1) continue;
                    std::vector<std::pair<internal::Vec3d, internal::Vec3d>> merged;
                    std::vector<bool> taken(cuts.size(), false);

                    for (size_t i = 0; i < cuts.size(); ++i) {
                        if (taken[i]) continue;
                        const internal::Vec3d& A = cuts[i].first;
                        const internal::Vec3d  AB = internal::sub(cuts[i].second, A);
                        const double           ab_len = internal::norm(AB);
                        if (ab_len < 1e-12) { taken[i] = true; continue; }
                        const internal::Vec3d  dir = internal::scale(AB, 1.0 / ab_len);

                        double t_lo = 0.0, t_hi = ab_len;
                        internal::Vec3d pt_lo = A;
                        internal::Vec3d pt_hi = cuts[i].second;
                        taken[i] = true;

                        for (size_t j = i + 1; j < cuts.size(); ++j) {
                            if (taken[j]) continue;
                            // Both endpoints of cuts[j] must lie on the line
                            // through A along dir.
                            auto on_line = [&](const internal::Vec3d& P) {
                                const internal::Vec3d AP = internal::sub(P, A);
                                const double t = internal::dot(AP, dir);
                                const internal::Vec3d proj = internal::add(
                                    A, internal::scale(dir, t));
                                const internal::Vec3d d = internal::sub(P, proj);
                                return internal::norm(d) < line_tol;
                            };
                            if (!on_line(cuts[j].first) || !on_line(cuts[j].second)) continue;

                            // Project j's endpoints to t-coordinates; extend
                            // [t_lo, t_hi] to cover them.
                            const double tj0 = internal::dot(internal::sub(cuts[j].first, A),  dir);
                            const double tj1 = internal::dot(internal::sub(cuts[j].second, A), dir);
                            if (tj0 < t_lo) { t_lo = tj0; pt_lo = cuts[j].first; }
                            if (tj0 > t_hi) { t_hi = tj0; pt_hi = cuts[j].first; }
                            if (tj1 < t_lo) { t_lo = tj1; pt_lo = cuts[j].second; }
                            if (tj1 > t_hi) { t_hi = tj1; pt_hi = cuts[j].second; }
                            taken[j] = true;
                        }
                        merged.emplace_back(pt_lo, pt_hi);
                    }
                    cuts.swap(merged);
                }
                size_t after_merge = 0;
                for (auto& v : tri_cuts) after_merge += v.size();
                SOUPT("after co-linear consolidation: %zu cut-records remain", after_merge);
            }

            // Split each filtered triangle by its accumulated cuts, then
            // re-classify each fragment via the voxel oracle. Build the
            // post-cut mesh.
            Mesh fragmented;
            fragmented.vertices.reserve(filtered.vertices.size() * 2);
            fragmented.faces.reserve(filtered.faces.size() * 2);

            // Position-based vertex dedup to avoid an explosion of
            // duplicate vertices from the cut endpoints (intersection
            // segments are shared across two triangles).
            struct Vec3Key {
                double x, y, z;
                bool operator==(const Vec3Key& o) const {
                    return x == o.x && y == o.y && z == o.z;
                }
            };
            struct Vec3KeyHash {
                std::size_t operator()(const Vec3Key& k) const {
                    auto h = std::hash<double>{}(k.x);
                    h ^= std::hash<double>{}(k.y) + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
                    h ^= std::hash<double>{}(k.z) + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
                    return h;
                }
            };
            std::unordered_map<Vec3Key, std::uint32_t, Vec3KeyHash> vmap;
            vmap.reserve(filtered.vertices.size() * 2);
            auto intern_vertex = [&](const internal::Vec3d& p) -> std::uint32_t {
                Vec3Key k{p[0], p[1], p[2]};
                auto it = vmap.find(k);
                if (it != vmap.end()) return it->second;
                const std::uint32_t id = static_cast<std::uint32_t>(fragmented.vertices.size());
                vmap.emplace(k, id);
                fragmented.vertices.push_back({p[0], p[1], p[2]});
                return id;
            };

            const double edge_tol = cell * 0.1;
            size_t total_pieces = 0;
            size_t total_kept = 0;

            for (std::uint32_t fi = 0; fi < filtered.faces.size(); ++fi) {
                const auto& f = filtered.faces[fi];
                const internal::Vec3d V0 = filtered.vertices[f[0]];
                const internal::Vec3d V1 = filtered.vertices[f[1]];
                const internal::Vec3d V2 = filtered.vertices[f[2]];

                std::vector<std::array<internal::Vec3d, 3>> pieces;
                if (tri_cuts[fi].empty()) {
                    pieces.push_back({V0, V1, V2});
                } else {
                    pieces = internal::split_triangle_by_cuts(
                        V0, V1, V2, tri_cuts[fi], edge_tol);
                }
                total_pieces += pieces.size();

                // Re-classify each piece via the voxel oracle. Keep
                // boundary candidates (one side Outside, the other
                // Inside), discard interior/exterior dust.
                for (const auto& p : pieces) {
                    const internal::Vec3d& A = p[0];
                    const internal::Vec3d& B = p[1];
                    const internal::Vec3d& C = p[2];
                    const std::array<double, 3> cen = {
                        (A[0]+B[0]+C[0]) / 3.0,
                        (A[1]+B[1]+C[1]) / 3.0,
                        (A[2]+B[2]+C[2]) / 3.0
                    };
                    const double e1[3] = { B[0]-A[0], B[1]-A[1], B[2]-A[2] };
                    const double e2[3] = { C[0]-A[0], C[1]-A[1], C[2]-A[2] };
                    double n[3] = {
                        e1[1]*e2[2] - e1[2]*e2[1],
                        e1[2]*e2[0] - e1[0]*e2[2],
                        e1[0]*e2[1] - e1[1]*e2[0]
                    };
                    const double nlen = std::sqrt(n[0]*n[0]+n[1]*n[1]+n[2]*n[2]);
                    if (nlen < 1e-30) continue;  // sliver from split
                    n[0] /= nlen; n[1] /= nlen; n[2] /= nlen;

                    using L = internal::VoxelGrid::Label;
                    const double sample_off = cell * 2.0;
                    const std::array<double, 3> p_pos = {
                        cen[0] + n[0]*sample_off, cen[1] + n[1]*sample_off, cen[2] + n[2]*sample_off
                    };
                    const std::array<double, 3> p_neg = {
                        cen[0] - n[0]*sample_off, cen[1] - n[1]*sample_off, cen[2] - n[2]*sample_off
                    };
                    const auto lp = grid.sample(p_pos);
                    const auto ln = grid.sample(p_neg);

                    bool keep = false;
                    if ((lp == L::Outside && ln == L::Inside) ||
                        (lp == L::Inside  && ln == L::Outside)) {
                        keep = true;
                    } else if (lp == L::Surface && ln == L::Surface) {
                        const double off2 = cell * 4.0;
                        const std::array<double, 3> q_pos = {
                            cen[0] + n[0]*off2, cen[1] + n[1]*off2, cen[2] + n[2]*off2
                        };
                        const std::array<double, 3> q_neg = {
                            cen[0] - n[0]*off2, cen[1] - n[1]*off2, cen[2] - n[2]*off2
                        };
                        const auto qlp = grid.sample(q_pos);
                        const auto qln = grid.sample(q_neg);
                        if (qlp != qln && qlp != L::Surface && qln != L::Surface) keep = true;
                    }
                    if (!keep) continue;

                    const std::uint32_t ia = intern_vertex(A);
                    const std::uint32_t ib = intern_vertex(B);
                    const std::uint32_t ic = intern_vertex(C);
                    if (ia == ib || ib == ic || ia == ic) continue;
                    fragmented.faces.push_back({ia, ib, ic});
                    ++total_kept;
                }
            }
            SOUPT("split: %zu input tris -> %zu pieces -> %zu kept",
                  filtered.faces.size(), total_pieces, total_kept);

            SOUPT("post-split+classify: %zu kept fragments (V=%zu)",
                  fragmented.faces.size(), fragmented.vertices.size());

            // Consolidate the post-cut mesh into a closed manifold via
            // welding (now at INPUT precision — fragments share endpoints
            // exactly via intern_vertex), orientation, and hole-fill.
            if (!fragmented.faces.empty()) {
                // Use a mild tolerance for the post-fragment weld; the
                // intersection endpoints are already unified by intern,
                // but original corner clusters in inputs like cube_soup
                // still need consolidation.
                auto wr = weld(fragmented, std::max(1e-9, g_diag * (1.0 / 30.0)));
                Mesh w = std::move(wr.mesh);

                auto or_ = orient_mesh(w);
                w = std::move(or_.mesh);

                auto hr = stages::fill_holes(w);
                w = std::move(hr.mesh);

                // Final NM cleanup + re-orient. Hole-fill across non-
                // edge-aligned fragment boundaries occasionally produces a
                // small number of duplicated faces around cube corners
                // that show up as NM edges. nm_edge's coherence-tiebreak
                // fin-removal usually resolves these.
                auto nme = stages::fix_non_manifold_edges(w);
                w = std::move(nme.mesh);
                if (nme.had_non_manifold_edges) {
                    auto or2 = orient_mesh(w);
                    w = std::move(or2.mesh);
                    auto hr2 = stages::fill_holes(w);
                    w = std::move(hr2.mesh);
                }

                auto vd = internal::compute_diagnostics(w);
                SOUPT("post weld+orient+holes+nm: V=%zu F=%zu bnd=%d nm=%d",
                      w.vertices.size(), w.faces.size(),
                      vd.open_boundary_edges, vd.non_manifold_edges);
                if (vd.open_boundary_edges == 0 && vd.non_manifold_edges == 0 &&
                    !w.faces.empty()) {
                    result.faces_after = static_cast<uint32_t>(w.faces.size());
                    result.success = true;
                    result.mesh = std::move(w);
                    return result;
                }
            }
            // Voxel-filter path didn't reach a clean manifold; fall through
            // to the legacy aggressive-weld + manifold::Manifold attempt
            // below using the ORIGINAL (un-filtered) mesh, since manifold
            // may still close a soup that the voxel filter mis-classified.
        }
    }

    // Compute bbox diagonal for soup-specific aggressive welding. Triangle
    // soup fixtures often have vertex positions perturbed by a few percent of
    // edge length (e.g., cube_soup has corners spread over ~0.04 of a unit
    // cube). The default scale-aware tolerance (bbox * 1e-7) is far too tight
    // to merge those clusters, so manifold sees 36 disconnected vertices and
    // bails out. We use a much larger soup tolerance — (bbox / 30) — which
    // collapses corner clusters but won't bridge distinct features in any
    // mesh worth a CSG repair.
    Mesh working = mesh;
    {
        double lo[3] = { std::numeric_limits<double>::max(),
                         std::numeric_limits<double>::max(),
                         std::numeric_limits<double>::max() };
        double hi[3] = { std::numeric_limits<double>::lowest(),
                         std::numeric_limits<double>::lowest(),
                         std::numeric_limits<double>::lowest() };
        for (const auto& v : working.vertices) {
            for (int k = 0; k < 3; ++k) {
                lo[k] = std::min(lo[k], v[k]);
                hi[k] = std::max(hi[k], v[k]);
            }
        }
        const double dx = hi[0] - lo[0];
        const double dy = hi[1] - lo[1];
        const double dz = hi[2] - lo[2];
        const double bbox_diag = std::sqrt(dx*dx + dy*dy + dz*dz);
        // Soup tolerance: ~3% of bbox diagonal, with a sane floor.
        const double soup_tol = std::max(1e-9, bbox_diag * (1.0 / 30.0));
        auto wr = weld(working, soup_tol);
        working = std::move(wr.mesh);
    }
    {
        auto or_ = orient_mesh(working);
        working = std::move(or_.mesh);
    }

    // Re-check: preprocessing may have produced a manifold mesh directly.
    auto diag2 = internal::compute_diagnostics(working);
    bool now_manifold = (diag2.open_boundary_edges == 0 && diag2.non_manifold_edges == 0);
    if (now_manifold) {
        result.faces_after = static_cast<uint32_t>(working.faces.size());
        result.success     = true;
        result.mesh        = std::move(working);
        return result;
    }

    // Re-compute bbox of `working` (post-weld) to set Manifold's tolerance.
    // The default (auto-computed) tolerance is too tight for triangle soup
    // inputs whose vertices are perturbed by a few percent of edge length.
    double w_lo[3] = { std::numeric_limits<double>::max(),
                       std::numeric_limits<double>::max(),
                       std::numeric_limits<double>::max() };
    double w_hi[3] = { std::numeric_limits<double>::lowest(),
                       std::numeric_limits<double>::lowest(),
                       std::numeric_limits<double>::lowest() };
    for (const auto& v : working.vertices) {
        for (int k = 0; k < 3; ++k) {
            w_lo[k] = std::min(w_lo[k], v[k]);
            w_hi[k] = std::max(w_hi[k], v[k]);
        }
    }
    const double w_dx = w_hi[0] - w_lo[0];
    const double w_dy = w_hi[1] - w_lo[1];
    const double w_dz = w_hi[2] - w_lo[2];
    const double w_diag = std::sqrt(w_dx*w_dx + w_dy*w_dy + w_dz*w_dz);

    // Build MeshGL from working mesh
    manifold::MeshGL mesh_gl;
    mesh_gl.numProp = 3;
    // Generous mesh-collapse tolerance for soup inputs. ~3% of bbox.
    mesh_gl.tolerance = static_cast<float>(std::max(1e-6, w_diag * 0.03));

    mesh_gl.vertProperties.reserve(working.vertices.size() * 3);
    for (const auto& v : working.vertices) {
        mesh_gl.vertProperties.push_back(static_cast<float>(v[0]));
        mesh_gl.vertProperties.push_back(static_cast<float>(v[1]));
        mesh_gl.vertProperties.push_back(static_cast<float>(v[2]));
    }

    mesh_gl.triVerts.reserve(working.faces.size() * 3);
    for (const auto& f : working.faces) {
        mesh_gl.triVerts.push_back(static_cast<uint32_t>(f[0]));
        mesh_gl.triVerts.push_back(static_cast<uint32_t>(f[1]));
        mesh_gl.triVerts.push_back(static_cast<uint32_t>(f[2]));
    }

    // Construct Manifold (soup import)
    manifold::Manifold m(mesh_gl);

    if (m.Status() != manifold::Manifold::Error::NoError || m.NumTri() == 0) {
        result.mesh = mesh;
        return result;
    }

    manifold::MeshGL64 out = m.GetMeshGL64();

    if (out.numProp < 3 || out.triVerts.empty() || out.vertProperties.empty()) {
        result.mesh = mesh;
        return result;
    }

    // Convert MeshGL64 back to Mesh (double precision, no float32 truncation)
    Mesh out_mesh;
    const uint64_t num_verts64 = static_cast<uint64_t>(out.vertProperties.size() / out.numProp);
    out_mesh.vertices.reserve(static_cast<size_t>(num_verts64));
    for (uint64_t i = 0; i < num_verts64; ++i) {
        uint64_t base = i * out.numProp;
        out_mesh.vertices.push_back({
            out.vertProperties[base + 0],
            out.vertProperties[base + 1],
            out.vertProperties[base + 2]
        });
    }

    const uint64_t num_tris64 = static_cast<uint64_t>(out.triVerts.size() / 3);
    out_mesh.faces.reserve(static_cast<size_t>(num_tris64));
    for (uint64_t i = 0; i < num_tris64; ++i) {
        uint64_t base = i * 3;
        out_mesh.faces.push_back({
            static_cast<uint32_t>(out.triVerts[base + 0]),
            static_cast<uint32_t>(out.triVerts[base + 1]),
            static_cast<uint32_t>(out.triVerts[base + 2])
        });
    }

    result.faces_after = static_cast<uint32_t>(out_mesh.faces.size());
    result.success     = true;
    result.mesh        = std::move(out_mesh);
    return result;
}

} // namespace meshseal::stages
