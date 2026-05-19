#include "shells.h"
#include "../internal/vec3.h"
#include <algorithm>
#include <cmath>
#include <queue>
#include <unordered_map>
#include <vector>

namespace meshseal::stages {

namespace {

static uint64_t undirected_edge_key(uint32_t a, uint32_t b) {
    if (a > b) std::swap(a, b);
    return (static_cast<uint64_t>(a) << 32) | static_cast<uint64_t>(b);
}

// Möller–Trumbore ray-triangle intersection along a configurable direction.
// Returns true if intersection found with t > epsilon along the +dir ray
// starting at origin. Previously hardcoded to {1,0,0}; now takes `dir` so
// the shells stage can do 3-ray majority voting for containment (originals'
// §6 prescription — robust against geometry tangent to the ray axis, which
// the single-ray version would mis-classify).
static bool ray_triangle_intersect(
    const internal::Vec3d& origin,
    const internal::Vec3d& dir,
    const internal::Vec3d& v0,
    const internal::Vec3d& v1,
    const internal::Vec3d& v2)
{
    constexpr double kParallelEps = 1e-10; // area-scale: cross-product magnitude
    constexpr double kTMinEps     = 1e-8;  // distance-scale: minimum valid t
    const internal::Vec3d edge1 = internal::sub(v1, v0);
    const internal::Vec3d edge2 = internal::sub(v2, v0);
    const internal::Vec3d h = internal::cross(dir, edge2);
    const double a = internal::dot(edge1, h);
    if (a > -kParallelEps && a < kParallelEps)
        return false; // ray parallel to triangle
    const double f = 1.0 / a;
    const internal::Vec3d s = internal::sub(origin, v0);
    const double u = f * internal::dot(s, h);
    if (u < 0.0 || u > 1.0)
        return false;
    const internal::Vec3d q = internal::cross(s, edge1);
    const double v = f * internal::dot(dir, q);
    if (v < 0.0 || (u + v) > 1.0)
        return false;
    const double t = f * internal::dot(edge2, q);
    return t > kTMinEps;
}

} // anonymous namespace

ShellResult analyze_shells(const Mesh& mesh, double vol_threshold) {
    ShellResult result;
    result.shells_dropped = 0;
    result.shells_kept    = 0;

    const uint32_t nfaces = static_cast<uint32_t>(mesh.faces.size());
    const uint32_t nverts = static_cast<uint32_t>(mesh.vertices.size());

    if (nfaces == 0) {
        result.mesh = mesh;
        return result;
    }

    // ----------------------------------------------------------------
    // Step 1: Build face adjacency via shared (undirected) edges, then BFS.
    // ----------------------------------------------------------------
    std::unordered_map<uint64_t, std::vector<uint32_t>> edge_to_faces;
    edge_to_faces.reserve(nfaces * 3);

    for (uint32_t fi = 0; fi < nfaces; ++fi) {
        const auto& f = mesh.faces[fi];
        for (int k = 0; k < 3; ++k) {
            uint32_t a = f[k], b = f[(k + 1) % 3];
            edge_to_faces[undirected_edge_key(a, b)].push_back(fi);
        }
    }

    // Build face neighbour list (only via manifold edges, i.e. exactly 2 faces).
    // For shells we allow boundary edges too; we connect faces sharing *any* edge.
    std::vector<std::vector<uint32_t>> face_adj(nfaces);
    for (const auto& kv : edge_to_faces) {
        const std::vector<uint32_t>& incident = kv.second;
        // Connect every pair of incident faces (handles manifold and non-manifold).
        for (size_t i = 0; i < incident.size(); ++i) {
            for (size_t j = i + 1; j < incident.size(); ++j) {
                face_adj[incident[i]].push_back(incident[j]);
                face_adj[incident[j]].push_back(incident[i]);
            }
        }
    }

    std::vector<int32_t> face_comp(nfaces, -1);
    int32_t num_comps = 0;
    std::queue<uint32_t> bfs;

    for (uint32_t seed = 0; seed < nfaces; ++seed) {
        if (face_comp[seed] != -1) continue;
        int32_t comp = num_comps++;
        face_comp[seed] = comp;
        bfs.push(seed);
        while (!bfs.empty()) {
            uint32_t cur = bfs.front();
            bfs.pop();
            for (uint32_t nb : face_adj[cur]) {
                if (face_comp[nb] == -1) {
                    face_comp[nb] = comp;
                    bfs.push(nb);
                }
            }
        }
    }

    // ----------------------------------------------------------------
    // Step 2: Compute signed volume per component.
    // ----------------------------------------------------------------
    std::vector<double> comp_volume(static_cast<size_t>(num_comps), 0.0);
    std::vector<uint32_t> comp_face_count(static_cast<size_t>(num_comps), 0u);

    for (uint32_t fi = 0; fi < nfaces; ++fi) {
        const auto& f = mesh.faces[fi];
        const internal::Vec3d& v0 = mesh.vertices[f[0]];
        const internal::Vec3d& v1 = mesh.vertices[f[1]];
        const internal::Vec3d& v2 = mesh.vertices[f[2]];
        int32_t comp = face_comp[fi];
        comp_volume[static_cast<size_t>(comp)] += internal::signed_tet_volume(v0, v1, v2);
        comp_face_count[static_cast<size_t>(comp)]++;
    }

    // ----------------------------------------------------------------
    // Step 3: Apply volume threshold.
    // ----------------------------------------------------------------
    double max_abs_vol = 0.0;
    for (int32_t c = 0; c < num_comps; ++c) {
        double av = std::abs(comp_volume[static_cast<size_t>(c)]);
        if (av > max_abs_vol) max_abs_vol = av;
    }

    std::vector<bool> comp_dropped(static_cast<size_t>(num_comps), false);
    // If all volumes are zero (e.g. intersecting soup with canceling windings),
    // skip the filter so downstream stages (soup_reconstruct) can handle it.
    if (max_abs_vol > 0.0) {
        for (int32_t c = 0; c < num_comps; ++c) {
            double av = std::abs(comp_volume[static_cast<size_t>(c)]);
            if (av == 0.0 || av < vol_threshold * max_abs_vol) {
                comp_dropped[static_cast<size_t>(c)] = true;
            }
        }
    }

    // ----------------------------------------------------------------
    // Step 4: Containment testing (best-effort, skip if > 100 shells).
    // ----------------------------------------------------------------
    // Collect kept component indices.
    std::vector<int32_t> kept_comps;
    for (int32_t c = 0; c < num_comps; ++c) {
        if (!comp_dropped[static_cast<size_t>(c)])
            kept_comps.push_back(c);
    }

    std::vector<bool> comp_contained(static_cast<size_t>(num_comps), false);

    if (kept_comps.size() > 1 && kept_comps.size() <= 100) {
        // For each kept component, compute a centroid of one of its faces.
        std::vector<internal::Vec3d> comp_centroid(static_cast<size_t>(num_comps), {0.0, 0.0, 0.0});
        std::vector<bool> centroid_set(static_cast<size_t>(num_comps), false);

        for (uint32_t fi = 0; fi < nfaces; ++fi) {
            int32_t c = face_comp[fi];
            if (comp_dropped[static_cast<size_t>(c)]) continue;
            if (centroid_set[static_cast<size_t>(c)]) continue;
            const auto& f = mesh.faces[fi];
            const internal::Vec3d& v0 = mesh.vertices[f[0]];
            const internal::Vec3d& v1 = mesh.vertices[f[1]];
            const internal::Vec3d& v2 = mesh.vertices[f[2]];
            comp_centroid[static_cast<size_t>(c)] = {
                (v0[0] + v1[0] + v2[0]) / 3.0,
                (v0[1] + v1[1] + v2[1]) / 3.0,
                (v0[2] + v1[2] + v2[2]) / 3.0
            };
            centroid_set[static_cast<size_t>(c)] = true;
        }

        // Pre-collect faces per kept component for ray casting against.
        // For each pair (outer, inner): cast ray from inner's centroid against outer's faces.
        for (size_t bi = 0; bi < kept_comps.size(); ++bi) {
            int32_t inner = kept_comps[bi];
            double inner_vol = std::abs(comp_volume[static_cast<size_t>(inner)]);

            for (size_t ai = 0; ai < kept_comps.size(); ++ai) {
                if (ai == bi) continue;
                int32_t outer = kept_comps[ai];
                double outer_vol = std::abs(comp_volume[static_cast<size_t>(outer)]);

                // Only test if outer is meaningfully larger than inner.
                if (outer_vol < inner_vol) continue;

                // Single-ray containment (originals' §6 prescribed 3-ray
                // majority vote; tried 2026-05-14 and changed sweep ±2 in
                // noisy ways without clear win — reverted, kept the
                // ray_triangle_intersect dir-parameter refactor for any
                // future re-attempt).
                const internal::Vec3d& origin = comp_centroid[static_cast<size_t>(inner)];
                const internal::Vec3d ray_dir = {1.0, 0.0, 0.0};
                int intersections = 0;
                for (uint32_t fi = 0; fi < nfaces; ++fi) {
                    if (face_comp[fi] != outer) continue;
                    const auto& f = mesh.faces[fi];
                    if (ray_triangle_intersect(
                            origin, ray_dir,
                            mesh.vertices[f[0]],
                            mesh.vertices[f[1]],
                            mesh.vertices[f[2]])) {
                        ++intersections;
                    }
                }
                if (intersections % 2 == 1) {
                    // Sign-aware refinement attempted 2026-05-13 (drop only
                    // on sign mismatch) — turned out backwards: same-sign
                    // nested = redundant duplicate (drop is fine), opposite-
                    // sign = legitimate hollow chamber (preserve). The fix
                    // would have dropped legitimate hollows. Reverted to
                    // unconditional drop; preservation gains from this stage
                    // come from the lowered vol_threshold (0.01 → 0.001).
                    comp_contained[static_cast<size_t>(inner)] = true;
                    comp_dropped[static_cast<size_t>(inner)]   = true;
                    break;
                }
            }
        }
    }

    // ----------------------------------------------------------------
    // Step 5: Build ShellInfo list.
    // ----------------------------------------------------------------
    for (int32_t c = 0; c < num_comps; ++c) {
        ShellInfo si;
        si.component_index = static_cast<uint32_t>(c);
        si.face_count      = comp_face_count[static_cast<size_t>(c)];
        si.volume          = comp_volume[static_cast<size_t>(c)];
        si.dropped         = comp_dropped[static_cast<size_t>(c)];
        si.is_contained    = comp_contained[static_cast<size_t>(c)];
        result.shells.push_back(si);

        if (si.dropped)
            ++result.shells_dropped;
        else
            ++result.shells_kept;
    }

    // ----------------------------------------------------------------
    // Step 6: Build output mesh — include only faces from kept components.
    // ----------------------------------------------------------------
    std::vector<std::array<uint32_t, 3>> new_faces;
    new_faces.reserve(nfaces);
    for (uint32_t fi = 0; fi < nfaces; ++fi) {
        int32_t c = face_comp[fi];
        if (!comp_dropped[static_cast<size_t>(c)])
            new_faces.push_back(mesh.faces[fi]);
    }

    // Compact vertices.
    std::vector<bool> used_vert(nverts, false);
    for (const auto& f : new_faces)
        for (int k = 0; k < 3; ++k)
            used_vert[f[k]] = true;

    std::vector<uint32_t> remap(nverts, UINT32_MAX);
    uint32_t new_idx = 0;
    std::vector<std::array<double, 3>> new_verts;
    new_verts.reserve(nverts);
    for (uint32_t vi = 0; vi < nverts; ++vi) {
        if (used_vert[vi]) {
            remap[vi] = new_idx++;
            new_verts.push_back(mesh.vertices[vi]);
        }
    }

    for (auto& f : new_faces)
        for (int k = 0; k < 3; ++k)
            f[k] = remap[f[k]];

    result.mesh.vertices = std::move(new_verts);
    result.mesh.faces    = std::move(new_faces);
    return result;
}

} // namespace meshseal::stages
