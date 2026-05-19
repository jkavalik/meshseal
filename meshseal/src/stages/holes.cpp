#include "holes.h"
#include "../internal/vec3.h"
#include "../internal/voxel_grid.h"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>
#include <list>
#include <memory>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace meshseal::stages {

using namespace meshseal::internal;

static uint64_t edge_key(uint32_t u, uint32_t v) {
    return (static_cast<uint64_t>(u) << 32) | static_cast<uint64_t>(v);
}

// Returns true if point p lies inside (or on the boundary of) triangle (a, b, c).
// Uses the sign consistency of cross products projected onto n.
static bool point_in_triangle_2d(
    const Vec3d& p,
    const Vec3d& a, const Vec3d& b, const Vec3d& c,
    const Vec3d& n)
{
    double d0 = dot(cross(sub(b, a), sub(p, a)), n);
    double d1 = dot(cross(sub(c, b), sub(p, b)), n);
    double d2 = dot(cross(sub(a, c), sub(p, c)), n);
    bool all_pos = (d0 >= 0.0) && (d1 >= 0.0) && (d2 >= 0.0);
    bool all_neg = (d0 <= 0.0) && (d1 <= 0.0) && (d2 <= 0.0);
    return all_pos || all_neg;
}

// ────────────────────────────────────────────────────────────────────────
// Liepa hole-fill — weighted minimum (max-dihedral, area) DP.
//
// Reference: Liepa, "Filling Holes in Meshes" (2003).
//
// For each candidate triangulation, the cost is a pair
//   (max dihedral over all internal/boundary edges, total area)
// compared lexicographically (max-dihedral dominant). This biases the DP
// toward fills that "follow the surface curvature" — adjacent fill
// triangles agree closely with each other and with the existing mesh
// triangles across the boundary edge they share. This is what avoids the
// failure mode the simple area-only Liepa hit on `doorman`: a valid
// triangulation whose triangles happen to cut through other parts of the
// model. The dihedral-aware cost rejects those when a smoother
// alternative exists.
//
// Triangle normals are computed in the FILL-ALIGNED direction (matching
// the existing mesh face's outward normal on the boundary edge), so
// dihedral score = 1 - dot(normal_a, normal_b) is 0 for coplanar / 2 for
// folded-back surfaces.
//
// Boundary-neighbor normals `bn[i]` give the outward normal of the
// existing mesh face on the boundary edge loop[i] → loop[(i+1)%n].
//
// Complexity: O(n^3) time, O(n^2) memory. Bounded by `max_n` (default 800)
// per loop; oversized loops fall through to the caller's ear-clipping
// fallback. Returns triangles in (loop[i], loop[k], loop[j]) order with
// i<k<j; caller flips winding to its convention.
struct LiepaCell {
    double max_dh    = 0.0;     // max dihedral score across edges; 0=flat
    double total_area = 0.0;
    int    K         = -1;      // split index (i < K < j)
    Vec3d  root_nrm  = {0,0,0}; // fill-aligned normal of triangle (i, K, j)
};

// Dihedral score: 1 - dot(n1, n2). 0 = coplanar same direction;
// 1 = perpendicular; 2 = anti-parallel (180° fold).
//
// NOTE: dihedral weighting is currently DISABLED (returns 0 for every pair).
// The DP infrastructure is in place — `bn` is computed, per-cell root
// normals are stored, the recurrence inspects both `dh_l` and `dh_r` — so
// re-enabling is a one-line change. The reason for the disable: with
// dihedral active, the (max-dihedral, area) lexicographic cost can pick a
// triangulation that physically crosses other parts of the model when
// adjacent geometry is highly curved. This was empirically caught by the
// `Pipeline: three chain-overlapping boxes` unit test, where weighted DP
// chose fills that the intersections-stage CSG then partially subtracted,
// dropping the unioned volume from 16 to 12. With dihedral disabled the
// DP picks pure min-area, identical to simple Liepa, and produces the
// same self-intersection count as the weighted version on the real-world
// fixtures we measured — so the dihedral term wasn't pulling its weight
// in the first place. A more discriminating weight (e.g. one that
// penalises triangles whose plane crosses neighbouring mesh tris,
// detectable via the voxel oracle we already have) is the right
// follow-up.
static inline double dihedral_score(const Vec3d& a, const Vec3d& b) {
    (void)a; (void)b;
    return 0.0;
}

// Voxel-oracle score for a candidate fill triangle. Sample the grid at
// centroid ± normal·(1.5·cell). If both samples have the same label, the
// triangle does NOT separate Inside from Outside — it either floats in the
// exterior (both Outside) or cuts through the interior (both Inside).
// Either way, it's a bad fill: penalty 1.0. Otherwise 0.0.
//
// Returns 0.0 if the grid is null (signal disabled) so call sites can be
// uniform. Centroid is the average of the three vertex positions in world
// coords; nrm is the fill-aligned unit normal as produced by
// tri_fill_normal().
//
// This is the modern replacement for dihedral weighting: instead of
// asking "is this triangle smooth with its neighbors?" (a local geometric
// signal) we ask "is this triangle in the right place to be a boundary?"
// (a global topological signal grounded in the input geometry).
// Multi-point: sample at 4 points across the triangle (centroid + 3 points
// each ~halfway toward a vertex). At each sample, probe ±normal*delta.
// Return the fraction of samples where both sides have the same Inside/Outside
// label (i.e. the triangle area locally fails to separate inside from outside).
// Score range: [0, 1]. Surface labels are skipped (neutral).
//
// Why multi-sample: Liepa's failure mode is fill triangles that graze interior
// geometry near a vertex while their centroid sits in clean space. Single-point
// centroid sampling misses these; sampling across the triangle catches them.
static inline double voxel_score(
    const Vec3d&            va,
    const Vec3d&            vb,
    const Vec3d&            vc,
    const Vec3d&            nrm,
    const internal::VoxelGrid* grid,
    double                  cell)
{
    if (!grid) return 0.0;
    const double mag = std::sqrt(dot(nrm, nrm));
    if (mag < 1e-30) return 0.0;
    const double inv = 1.0 / mag;
    const double delta = cell * 2.0;
    const Vec3d un = {nrm[0]*inv, nrm[1]*inv, nrm[2]*inv};

    // 4 sample positions on the triangle: centroid + 3 near-vertex points
    // (each is the centroid of {vertex, centroid, centroid} = 1/3 toward vertex).
    const Vec3d cen = {(va[0]+vb[0]+vc[0])/3.0,
                       (va[1]+vb[1]+vc[1])/3.0,
                       (va[2]+vb[2]+vc[2])/3.0};
    auto mid = [&](const Vec3d& v){
        return Vec3d{(v[0]+2.0*cen[0])/3.0,
                     (v[1]+2.0*cen[1])/3.0,
                     (v[2]+2.0*cen[2])/3.0};
    };
    const Vec3d samples[4] = {cen, mid(va), mid(vb), mid(vc)};

    using L = internal::VoxelGrid::Label;
    int bad = 0;
    int counted = 0;
    for (const auto& s : samples) {
        const std::array<double, 3> p_plus  = {s[0]+un[0]*delta, s[1]+un[1]*delta, s[2]+un[2]*delta};
        const std::array<double, 3> p_minus = {s[0]-un[0]*delta, s[1]-un[1]*delta, s[2]-un[2]*delta};
        const L a = grid->sample(p_plus);
        const L b = grid->sample(p_minus);
        if (a == L::Surface || b == L::Surface) continue;
        ++counted;
        if (a == b) ++bad;
    }
    if (counted == 0) return 0.0;
    return static_cast<double>(bad) / static_cast<double>(counted);
}

// Fill-aligned normal of triangle (loop[a], loop[k], loop[b]).
// Equal to cross((loop[b]-loop[a]), (loop[k]-loop[a])) normalised — note
// the swapped order; this aligns with the FILL side of the boundary (the
// triangle is emitted with reversed winding (loop[i], loop[j], loop[k]) by
// the caller, and this normal matches that orientation).
static Vec3d tri_fill_normal(
    int a, int k, int b,
    const std::vector<uint32_t>& loop,
    const std::vector<Vec3d>&    verts)
{
    const Vec3d& va = verts[loop[a]];
    const Vec3d& vk = verts[loop[k]];
    const Vec3d& vb = verts[loop[b]];
    const Vec3d nrm = cross(sub(vb, va), sub(vk, va));
    const double m = std::sqrt(dot(nrm, nrm));
    if (m < 1e-30) return {0.0, 0.0, 0.0};
    return {nrm[0] / m, nrm[1] / m, nrm[2] / m};
}

static std::vector<std::array<uint32_t, 3>> liepa_dp_triangulate(
    const std::vector<uint32_t>& loop,
    const std::vector<Vec3d>&    verts,
    const std::vector<Vec3d>&    bn,   // bn[i] = outward normal of mesh face on edge loop[i]→loop[(i+1)%n]
    const internal::VoxelGrid*   grid  = nullptr,
    double                       cell  = 0.0,
    int                          max_n = 800)
{
    const int n = static_cast<int>(loop.size());
    if (n < 3) return {};
    if (n > max_n) return {};
    if (n == 3) {
        return {{loop[0], loop[1], loop[2]}};
    }

    auto tri_area = [&](int a, int b, int c) -> double {
        const Vec3d& va = verts[loop[a]];
        const Vec3d& vb = verts[loop[b]];
        const Vec3d& vc = verts[loop[c]];
        const Vec3d eab = sub(vb, va);
        const Vec3d ebc = sub(vc, vb);
        const Vec3d eca = sub(va, vc);
        const Vec3d nrm = cross(eab, sub(vc, va));
        const double a2 = dot(nrm, nrm);
        // Heavy penalty for degenerate (near-zero-area) triangles —
        // otherwise the min-area DP picks them eagerly and they then get
        // filtered by remove_degenerate_faces, re-opening the boundary.
        if (a2 < 1e-40) return 1e20;
        const double area = 0.5 * std::sqrt(a2);
        // Aspect-ratio penalty (2026-05-14): on near-collinear or thin loops
        // Liepa's area-only weight happily picks sliver triangles (longest
        // edge >> shortest altitude) that downstream stages then drop or
        // mis-handle. Penalise slivers in the DP weight. Aspect quality
        // metric: `4 * sqrt(3) * area / (sum_of_edge_lengths_squared)`,
        // which equals 1 for an equilateral triangle and approaches 0 for
        // slivers (= Pillinger's "quality" measure). Penalty multiplier:
        // 1 + alpha * (1 - quality), with alpha small enough not to override
        // legitimate large-area choices. alpha = 0.5 → equilateral pays 1×,
        // a 1:10 sliver pays ~1.45× — biases toward regular triangulations
        // when areas are otherwise comparable.
        const double l2_ab = dot(eab, eab);
        const double l2_bc = dot(ebc, ebc);
        const double l2_ca = dot(eca, eca);
        const double sum_l2 = l2_ab + l2_bc + l2_ca;
        const double quality = (sum_l2 > 0.0)
            ? (4.0 * 1.7320508075688772 * area) / sum_l2
            : 0.0;
        const double aspect_penalty = 1.0 + 0.15 * (1.0 - std::min(1.0, std::max(0.0, quality)));
        return area * aspect_penalty;
    };

    std::vector<std::vector<LiepaCell>> dp(n, std::vector<LiepaCell>(n));

    // Lexicographic comparator: prefer smaller max_dh; ties broken by area.
    auto better = [](double a_dh, double a_area, double b_dh, double b_area) {
        if (std::abs(a_dh - b_dh) > 1e-12) return a_dh < b_dh;
        return a_area < b_area;
    };

    // Base case: chain length 2 (single triangle loop[i], loop[i+1], loop[i+2]).
    // Edges (i, i+1) and (i+1, i+2) are boundary edges with known mesh
    // neighbors bn[i] and bn[i+1]. Edge (i, i+2) is the closing diagonal,
    // handled by the parent.
    for (int i = 0; i + 2 < n; ++i) {
        const int j = i + 2;
        const int k = i + 1;
        const Vec3d nrm = tri_fill_normal(i, k, j, loop, verts);
        const double dh_left  = dihedral_score(nrm, bn[i]);
        const double dh_right = dihedral_score(nrm, bn[k]);
        const double vs = voxel_score(verts[loop[i]], verts[loop[k]], verts[loop[j]],
                                      nrm, grid, cell);
        dp[i][j] = {std::max({dh_left, dh_right, vs}), tri_area(i, k, j), k, nrm};
    }

    // Fill DP for chain lengths 3..n-1.
    for (int len = 3; len < n; ++len) {
        for (int i = 0; i + len < n; ++i) {
            const int j = i + len;
            double best_dh   = std::numeric_limits<double>::infinity();
            double best_area = std::numeric_limits<double>::infinity();
            int    best_k    = -1;
            Vec3d  best_nrm  = {0, 0, 0};
            for (int k = i + 1; k < j; ++k) {
                const Vec3d nrm = tri_fill_normal(i, k, j, loop, verts);
                // Left neighbor across edge (i, k):
                //   boundary edge if k == i+1 → mesh face normal bn[i]
                //   otherwise → root triangle normal of sub-DP dp[i][k]
                const Vec3d left_nbr =
                    (k == i + 1) ? bn[i] : dp[i][k].root_nrm;
                // Right neighbor across edge (k, j):
                const Vec3d right_nbr =
                    (j == k + 1) ? bn[k] : dp[k][j].root_nrm;

                double dh_l = dihedral_score(nrm, left_nbr);
                double dh_r = dihedral_score(nrm, right_nbr);
                const double vs = voxel_score(verts[loop[i]], verts[loop[k]], verts[loop[j]],
                                              nrm, grid, cell);

                double max_dh = std::max({dp[i][k].max_dh, dp[k][j].max_dh,
                                          dh_l, dh_r, vs});

                // Top-level closing edge (0, n-1) is the wrap-around polygon
                // edge from loop[n-1] to loop[0], so its neighbor is bn[n-1].
                if (i == 0 && j == n - 1) {
                    max_dh = std::max(max_dh, dihedral_score(nrm, bn[n - 1]));
                }

                double total_area = dp[i][k].total_area + dp[k][j].total_area +
                                    tri_area(i, k, j);

                if (better(max_dh, total_area, best_dh, best_area)) {
                    best_dh   = max_dh;
                    best_area = total_area;
                    best_k    = k;
                    best_nrm  = nrm;
                }
            }
            dp[i][j] = {best_dh, best_area, best_k, best_nrm};
        }
    }

    // Reconstruct triangles via iterative DFS over the K table.
    std::vector<std::array<uint32_t, 3>> tris;
    tris.reserve(static_cast<size_t>(n - 2));
    std::vector<std::pair<int, int>> stack;
    stack.emplace_back(0, n - 1);
    while (!stack.empty()) {
        const auto [i, j] = stack.back();
        stack.pop_back();
        if (j - i < 2) continue;
        const int k = dp[i][j].K;
        if (k < 0) continue;
        tris.push_back({loop[i], loop[k], loop[j]});
        stack.emplace_back(i, k);
        stack.emplace_back(k, j);
    }
    return tris;
}

HoleFillResult fill_holes(const Mesh& input, uint32_t max_fan_size) {
    HoleFillResult result;
    result.mesh         = input;
    result.holes_filled = 0;
    result.faces_added  = 0;

    // Use the original vertex positions for all lookups.
    // New centroid vertices are appended to result.mesh.vertices but never
    // need to be looked up by index during hole computation.
    const auto& verts       = input.vertices;
    const auto& input_faces = input.faces;

    // ---- Step 1: Count directed half-edges and derive boundary half-edges ----
    // For each undirected edge, the boundary contribution in direction u->v is
    //   max(0, count(u->v) - count(v->u))
    // i.e. the number of u->v half-edges that have no v->u twin. This handles
    // figure-8 boundaries (one vertex with multiple outgoing boundary edges)
    // and non-manifold edges with unbalanced winding correctly.
    std::unordered_map<uint64_t, uint32_t> directed_edge_count;
    directed_edge_count.reserve(input_faces.size() * 6);

    // Directed-edge → face-index map. Used to find the existing mesh face
    // on the OTHER side of each boundary edge, whose outward normal becomes
    // the boundary-neighbor normal for weighted Liepa.
    std::unordered_map<uint64_t, uint32_t> directed_to_face;
    directed_to_face.reserve(input_faces.size() * 3);

    for (uint32_t fi = 0; fi < input_faces.size(); ++fi) {
        const auto& f = input_faces[fi];
        directed_edge_count[edge_key(f[0], f[1])]++;
        directed_edge_count[edge_key(f[1], f[2])]++;
        directed_edge_count[edge_key(f[2], f[0])]++;
        directed_to_face[edge_key(f[0], f[1])] = fi;
        directed_to_face[edge_key(f[1], f[2])] = fi;
        directed_to_face[edge_key(f[2], f[0])] = fi;
    }

    // boundary_outgoing[u] = vector of v's such that there are unpaired u->v
    // half-edges, with multiplicity. Tracing each loop pops one outgoing edge
    // per visit, so vertices on multiple loops (figure-8 pinches) are handled.
    std::unordered_map<uint32_t, std::vector<uint32_t>> boundary_outgoing;
    for (const auto& kv : directed_edge_count) {
        uint64_t k = kv.first;
        uint32_t u = static_cast<uint32_t>(k >> 32);
        uint32_t v = static_cast<uint32_t>(k & 0xFFFFFFFFULL);
        const uint32_t cnt_uv = kv.second;
        auto it = directed_edge_count.find(edge_key(v, u));
        const uint32_t cnt_vu = (it == directed_edge_count.end()) ? 0u : it->second;
        if (cnt_uv > cnt_vu) {
            const uint32_t excess = cnt_uv - cnt_vu;
            auto& vs = boundary_outgoing[u];
            for (uint32_t i = 0; i < excess; ++i) vs.push_back(v);
        }
    }

    if (boundary_outgoing.empty()) {
        return result; // already closed (or empty)
    }

    // ---- Step 2: Trace closed boundary loops by consuming outgoing edges ----
    auto pop_outgoing = [&](uint32_t u) -> uint32_t {
        auto it = boundary_outgoing.find(u);
        if (it == boundary_outgoing.end() || it->second.empty()) return UINT32_MAX;
        uint32_t v = it->second.back();
        it->second.pop_back();
        return v;
    };
    auto any_remaining_start = [&]() -> uint32_t {
        for (auto& kv : boundary_outgoing) {
            if (!kv.second.empty()) return kv.first;
        }
        return UINT32_MAX;
    };

    std::vector<std::vector<uint32_t>> loops;
    while (true) {
        const uint32_t start = any_remaining_start();
        if (start == UINT32_MAX) break;

        std::vector<uint32_t> loop;
        uint32_t cur = start;
        bool closed = false;
        // Safety bound: each step consumes one boundary edge. Total boundary
        // edges <= 3 * |F|, so a single loop trace can't exceed that.
        const size_t max_steps = input_faces.size() * 3 + 1;
        for (size_t step = 0; step < max_steps; ++step) {
            loop.push_back(cur);
            uint32_t nxt = pop_outgoing(cur);
            if (nxt == UINT32_MAX) break;             // open chain — abandoned
            if (nxt == start) { closed = true; break; }
            cur = nxt;
        }
        if (closed && loop.size() >= 3) {
            loops.push_back(std::move(loop));
        }
        // Else: open chain; its consumed edges are already removed, so the
        // outer while loop continues with whatever remains.
    }

    // ---- Step 2b: Build voxel oracle for fill-placement scoring ----
    //
    // Empirically disabled (see liepa call site below). Flip this flag to
    // build the grid; flip the call site's nullptr to grid.get() to use it.
    constexpr bool kBuildVoxelOracle = false;
    //
    // Rasterize every input face, then seal each boundary loop with a
    // coarse centroid-fan rasterization (purely temporary, voxel-only —
    // these fans are NOT emitted as fill triangles). The seal closes the
    // holes in voxel space so flood_fill_exterior cannot leak into the
    // interior. After flood-fill, sample ± normal at each candidate fill
    // triangle's centroid to detect "fin through interior" or "floats in
    // air" placements; see voxel_score().
    //
    // Resolution 64 along the longest axis is a balance: coarse enough to
    // tolerate float wobble and small geometric noise without misclassifying,
    // fine enough that thin features (~bbox/64 thick) still register.
    std::unique_ptr<internal::VoxelGrid> grid;
    double grid_cell = 0.0;
    if (kBuildVoxelOracle) {
        if (!verts.empty()) {
            std::array<double, 3> lo = verts[0];
            std::array<double, 3> hi = verts[0];
            for (const auto& v : verts) {
                for (int k = 0; k < 3; ++k) {
                    lo[k] = std::min(lo[k], v[k]);
                    hi[k] = std::max(hi[k], v[k]);
                }
            }
            const double ext = std::max({hi[0]-lo[0], hi[1]-lo[1], hi[2]-lo[2]});
            if (ext > 0.0) {
                grid = std::make_unique<internal::VoxelGrid>(lo, hi, /*target=*/64);
                grid_cell = grid->cell_size();
                // Rasterize input faces.
                for (const auto& f : input_faces) {
                    grid->rasterize(verts[f[0]], verts[f[1]], verts[f[2]]);
                }
                // Rasterize centroid-fan triangulation of each loop (seal).
                for (const auto& loop : loops) {
                    const uint32_t n = static_cast<uint32_t>(loop.size());
                    Vec3d c = {0.0, 0.0, 0.0};
                    for (uint32_t idx : loop) c = add(c, verts[idx]);
                    c = scale(c, 1.0 / static_cast<double>(n));
                    const std::array<double, 3> cen = {c[0], c[1], c[2]};
                    for (uint32_t i = 0; i < n; ++i) {
                        const uint32_t ni = (i + 1) % n;
                        grid->rasterize(verts[loop[i]], verts[loop[ni]], cen);
                    }
                }
                grid->flood_fill_exterior();
            }
        }
    }

    // ---- Step 3: Fill each loop ----
    std::vector<std::array<uint32_t, 3>> new_faces;

    for (const auto& loop : loops) {
        const uint32_t n = static_cast<uint32_t>(loop.size());

        // Centroid of boundary loop vertices
        Vec3d centroid = {0.0, 0.0, 0.0};
        for (uint32_t idx : loop) {
            centroid = add(centroid, verts[idx]);
        }
        centroid = scale(centroid, 1.0 / static_cast<double>(n));

        // Loop normal = sum of fan cross-products from centroid.
        // Fill triangles must have normals OPPOSITE to this direction.
        // |loop_normal_sum| = 2 * (planar projected area of the loop).
        Vec3d loop_normal_sum = {0.0, 0.0, 0.0};
        for (uint32_t i = 0; i < n; i++) {
            Vec3d vi = sub(verts[loop[i]],           centroid);
            Vec3d vn = sub(verts[loop[(i + 1) % n]], centroid);
            loop_normal_sum = add(loop_normal_sum, cross(vi, vn));
        }
        Vec3d loop_normal = normalize(loop_normal_sum);

        // Guard against collinear / near-collinear loops: when all loop
        // vertices lie on a single line, the centroid lies on the same line,
        // so every fan triangle (centroid, v_i, v_{i+1}) is three-collinear
        // and has zero area. The fill would be silently dropped by
        // remove_degenerate_faces and the boundary would re-open. Detect this
        // early via the loop's bbox and area magnitude, and skip the fill;
        // the component will be preserved in the outer pipeline's pile.
        {
            double bbox_lo[3] = { verts[loop[0]][0], verts[loop[0]][1], verts[loop[0]][2] };
            double bbox_hi[3] = { bbox_lo[0], bbox_lo[1], bbox_lo[2] };
            for (uint32_t idx : loop) {
                for (int k = 0; k < 3; ++k) {
                    bbox_lo[k] = std::min(bbox_lo[k], verts[idx][k]);
                    bbox_hi[k] = std::max(bbox_hi[k], verts[idx][k]);
                }
            }
            const double dx = bbox_hi[0] - bbox_lo[0];
            const double dy = bbox_hi[1] - bbox_lo[1];
            const double dz = bbox_hi[2] - bbox_lo[2];
            const double scale_sq = dx * dx + dy * dy + dz * dz;
            const double area_mag = std::sqrt(dot(loop_normal_sum, loop_normal_sum));
            // For a non-degenerate planar loop with bbox diagonal d, area_mag
            // scales with d^2. Collinear loops have area_mag == 0; very thin
            // (slit-like) loops have area_mag ~ d * width. Reject anything
            // below 1e-10 of the squared scale (well below the 1e-14 area
            // threshold used by remove_degenerate_faces).
            if (scale_sq > 0.0 && area_mag < 1e-10 * scale_sq) {
                continue; // skip degenerate loop; boundary stays open
            }
        }

        const uint32_t faces_before = static_cast<uint32_t>(new_faces.size());

        // Boundary-neighbor normals: for each i, bn[i] is the outward normal
        // of the mesh face containing the directed boundary edge
        // loop[i] → loop[(i+1)%n]. Used by weighted Liepa to compute
        // dihedral with the existing mesh across boundary edges and to bias
        // the fill toward smoothness.
        std::vector<Vec3d> bn(n, Vec3d{0.0, 0.0, 0.0});
        for (uint32_t i = 0; i < n; ++i) {
            const uint32_t ni = (i + 1) % n;
            auto it = directed_to_face.find(edge_key(loop[i], loop[ni]));
            if (it == directed_to_face.end()) continue;
            const auto& f = input_faces[it->second];
            const Vec3d nrm = cross(sub(verts[f[1]], verts[f[0]]),
                                    sub(verts[f[2]], verts[f[0]]));
            const double m = std::sqrt(dot(nrm, nrm));
            if (m < 1e-30) continue;
            bn[i] = {nrm[0] / m, nrm[1] / m, nrm[2] / m};
        }

        // Primary path: Liepa weighted (max-dihedral, area) DP. No centroid
        // → no off-surface vertex → no fan-into-body self-intersections.
        // Up to ~800-vertex loops (O(n^3) DP); larger falls through to
        // ear-clipping.  The `max_fan_size` parameter is preserved for ABI
        // compatibility but no longer gates fan-fill (which is removed
        // entirely).
        (void)max_fan_size;
        // Voxel-aware weight DISABLED: empirically, the centroid+sub-centroid
        // multi-point signal didn't reduce SI on doorman/planar_mesh/black_vase
        // (numbers stayed within ±1 % of the disabled baseline; oracle is 4-5×
        // better). The remaining SI on real-world fixtures doesn't come from
        // candidates the DP can discriminate between — it's from the
        // intersections-stage CSG output, fin-removal cascades, or coplanar
        // overlaps. Infrastructure (grid build, voxel_score) preserved for a
        // smarter follow-up signal; flip the second arg back to `grid.get()`
        // to re-enable. See git history / CLAUDE.md for measurements.
        (void)grid; (void)grid_cell;
        auto liepa = liepa_dp_triangulate(loop, verts, bn, /*grid=*/nullptr, /*cell=*/0.0);
        bool filled = false;

        if (!liepa.empty()) {
            // Liepa returns triangles in (i, k, j) index order. For the
            // boundary trace direction we use, those triangles' natural
            // cross products align with +loop_normal_sum, so we need to
            // flip the winding to land on -loop_normal_sum (the mesh
            // outward direction across the boundary). Apply uniformly per
            // loop — Liepa is consistent within one loop.
            bool need_flip = true;
            {
                const auto& t = liepa.front();
                const Vec3d& va = verts[t[0]];
                const Vec3d& vb = verts[t[1]];
                const Vec3d& vc = verts[t[2]];
                const Vec3d nrm = cross(sub(vb, va), sub(vc, va));
                if (dot(nrm, loop_normal_sum) < 0.0) need_flip = false;
            }
            for (auto& t : liepa) {
                if (need_flip) std::swap(t[1], t[2]);
                new_faces.push_back(t);
            }
            filled = true;
            (void)centroid;
        }

        if (!filled) {
            // Fallback: ear-clipping for very large loops Liepa skipped.
            // Same algorithm as before this change.
            std::list<uint32_t> remaining(loop.begin(), loop.end());

            while (remaining.size() > 3) {
                bool ear_found = false;

                for (auto it = remaining.begin(); it != remaining.end(); ++it) {
                    auto prev_it = (it == remaining.begin())
                                       ? std::prev(remaining.end())
                                       : std::prev(it);
                    auto next_it = std::next(it);
                    if (next_it == remaining.end()) next_it = remaining.begin();

                    uint32_t ip  = *prev_it;
                    uint32_t iv  = *it;
                    uint32_t in_ = *next_it;

                    const Vec3d& vp = verts[ip];
                    const Vec3d& vv = verts[iv];
                    const Vec3d& vn = verts[in_];

                    Vec3d ear_cross = cross(sub(vv, vp), sub(vn, vp));
                    if (dot(ear_cross, loop_normal_sum) <= 0.0) continue;

                    bool has_inside = false;
                    for (auto ot = remaining.begin(); ot != remaining.end(); ++ot) {
                        if (*ot == ip || *ot == iv || *ot == in_) continue;
                        if (point_in_triangle_2d(verts[*ot], vp, vv, vn, loop_normal)) {
                            has_inside = true;
                            break;
                        }
                    }
                    if (has_inside) continue;

                    new_faces.push_back({in_, iv, ip});
                    remaining.erase(it);
                    ear_found = true;
                    break;
                }

                if (!ear_found) break;
            }

            if (remaining.size() == 3) {
                auto it = remaining.begin();
                uint32_t v0 = *it++;
                uint32_t v1 = *it++;
                uint32_t v2 = *it;
                new_faces.push_back({v2, v1, v0});
            }
        }

        const uint32_t faces_added_this =
            static_cast<uint32_t>(new_faces.size()) - faces_before;
        if (faces_added_this > 0) {
            result.holes_filled++;
            result.faces_added += faces_added_this;
        }
    }

    // Append all new faces to the mesh
    for (const auto& f : new_faces) {
        result.mesh.faces.push_back(f);
    }

    return result;
}

} // namespace meshseal::stages
