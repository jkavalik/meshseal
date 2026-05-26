#include "nm_edge.h"
#include "../internal/vec3.h"
#include <algorithm>
#include <cmath>
#include <numeric>
#include <queue>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace meshseal::stages {

namespace {

static uint64_t undirected_edge_key(uint32_t a, uint32_t b) {
    if (a > b) std::swap(a, b);
    return (static_cast<uint64_t>(a) << 32) | static_cast<uint64_t>(b);
}

} // anonymous namespace

NmEdgeResult fix_non_manifold_edges(const Mesh& mesh) {
    NmEdgeResult result;
    result.faces_removed = 0;
    result.edges_fixed = 0;
    result.had_non_manifold_edges = false;

    const uint32_t nfaces = static_cast<uint32_t>(mesh.faces.size());
    const uint32_t nverts = static_cast<uint32_t>(mesh.vertices.size());

    if (nfaces == 0) {
        result.mesh = mesh;
        return result;
    }

    // Step 1: Build edge_to_faces (undirected edge -> list of incident face indices)
    std::unordered_map<uint64_t, std::vector<uint32_t>> edge_to_faces;
    edge_to_faces.reserve(nfaces * 3);

    for (uint32_t fi = 0; fi < nfaces; ++fi) {
        const auto& f = mesh.faces[fi];
        for (int k = 0; k < 3; ++k) {
            uint32_t a = f[k], b = f[(k + 1) % 3];
            uint64_t key = undirected_edge_key(a, b);
            edge_to_faces[key].push_back(fi);
        }
    }

    // Step 2: Find all non-manifold edges (> 2 incident faces).
    // NOTE: iteration order is unordered_map-bucket-defined. This is a
    // documented non-determinism risk (sweep results CAN differ across
    // stdlib versions), but sorting here REGRESSES black_vase from CLEAN
    // to nm=1: the downstream nm_carve_refill rescue depends on the
    // specific pair choices the bucket order produces. Until the
    // downstream is decoupled from this order, leave as-is.
    std::vector<uint64_t> nm_edges;
    for (const auto& kv : edge_to_faces) {
        if (kv.second.size() > 2) {
            nm_edges.push_back(kv.first);
        }
    }

    if (nm_edges.empty()) {
        result.mesh = mesh;
        return result;
    }

    result.had_non_manifold_edges = true;

    // Step 3: Build connected components using only manifold edges (exactly 2 incident faces)
    // Faces adjacent only through non-manifold or boundary edges are isolated components.
    std::vector<std::vector<uint32_t>> face_neighbors(nfaces);
    for (const auto& kv : edge_to_faces) {
        if (kv.second.size() == 2) {
            uint32_t fa = kv.second[0];
            uint32_t fb = kv.second[1];
            face_neighbors[fa].push_back(fb);
            face_neighbors[fb].push_back(fa);
        }
    }

    std::vector<int32_t> face_component(nfaces, -1);
    std::vector<uint32_t> component_size;
    std::queue<uint32_t> bfs_queue;
    int32_t num_components = 0;

    for (uint32_t fi = 0; fi < nfaces; ++fi) {
        if (face_component[fi] != -1) continue;
        int32_t comp = num_components++;
        component_size.push_back(0);
        face_component[fi] = comp;
        bfs_queue.push(fi);
        while (!bfs_queue.empty()) {
            uint32_t cur = bfs_queue.front();
            bfs_queue.pop();
            component_size[comp]++;
            for (uint32_t nb : face_neighbors[cur]) {
                if (face_component[nb] == -1) {
                    face_component[nb] = comp;
                    bfs_queue.push(nb);
                }
            }
        }
    }

    // Step 3.5: Pre-compute unit face normals for the geometric tie-break.
    auto face_normal = [&](uint32_t fi) -> internal::Vec3d {
        const auto& f = mesh.faces[fi];
        const auto& v0 = mesh.vertices[f[0]];
        const auto& v1 = mesh.vertices[f[1]];
        const auto& v2 = mesh.vertices[f[2]];
        internal::Vec3d e1 = internal::sub(v1, v0);
        internal::Vec3d e2 = internal::sub(v2, v0);
        internal::Vec3d n = internal::cross(e1, e2);
        const double mag = internal::norm(n);
        if (mag < 1e-30) return {0.0, 0.0, 0.0};
        return internal::scale(n, 1.0 / mag);
    };

    // Step 4: For each non-manifold edge, keep one face per winding group;
    // tie-break by *normal coherence* with the face's non-NM neighbors so that
    // the face most consistent with its surrounding surface is preserved.
    // Index-based fallback only when no neighbor information is available.
    std::vector<bool> remove_face(nfaces, false);

    auto coherence_score = [&](uint32_t fi, uint64_t skip_edge_key) -> double {
        // Sum of dot(n_fi, n_neighbor) over neighbors via this face's two
        // OTHER edges (not skip_edge_key). Higher = better fit with surface.
        // Returns -DBL_MAX if no manifold neighbors are available.
        const auto& f = mesh.faces[fi];
        const internal::Vec3d nf = face_normal(fi);
        double sum = 0.0;
        int cnt = 0;
        for (int k = 0; k < 3; ++k) {
            uint32_t a = f[k], b = f[(k + 1) % 3];
            uint64_t ek = undirected_edge_key(a, b);
            if (ek == skip_edge_key) continue;
            auto it = edge_to_faces.find(ek);
            if (it == edge_to_faces.end()) continue;
            const auto& nbs = it->second;
            // Only score against manifold-adjacent neighbors (valence 2).
            if (nbs.size() != 2) continue;
            uint32_t other = (nbs[0] == fi) ? nbs[1] : nbs[0];
            if (remove_face[other]) continue;
            const internal::Vec3d nn = face_normal(other);
            sum += internal::dot(nf, nn);
            ++cnt;
        }
        if (cnt == 0) return -1e300;
        return sum / static_cast<double>(cnt);
    };

    for (size_t ei = 0; ei < nm_edges.size(); ++ei) {
        uint64_t ek = nm_edges[ei];
        const std::vector<uint32_t>& faces_for_edge = edge_to_faces[ek];

        // Decode endpoints: undirected_edge_key stores min << 32 | max
        uint32_t u = static_cast<uint32_t>(ek >> 32);
        uint32_t v = static_cast<uint32_t>(ek & 0xFFFFFFFF);

        // Separate faces into group A (directed u->v) and group B (directed v->u)
        std::vector<uint32_t> groupA, groupB;
        for (size_t i = 0; i < faces_for_edge.size(); ++i) {
            uint32_t fi = faces_for_edge[i];
            if (remove_face[fi]) continue; // already marked by a prior nm_edge
            const auto& f = mesh.faces[fi];
            bool found_uv = false;
            for (int k = 0; k < 3; ++k) {
                if (f[k] == u && f[(k + 1) % 3] == v) {
                    found_uv = true;
                    break;
                }
            }
            if (found_uv)
                groupA.push_back(fi);
            else
                groupB.push_back(fi);
        }

        // For each group with > 1 face, keep the most surface-coherent face
        // (highest average dot product with its non-NM-edge neighbors). Falls
        // back to "largest component, smallest index" when no manifold
        // neighbor info is available (e.g. all faces are isolated).
        auto mark_excess = [&](std::vector<uint32_t>& group) {
            if (group.size() <= 1) return;
            uint32_t best_fi = group[0];
            double best_coh = coherence_score(best_fi, ek);
            uint32_t best_cs = component_size[static_cast<size_t>(face_component[best_fi])];
            for (size_t i = 1; i < group.size(); ++i) {
                uint32_t fi = group[i];
                double coh = coherence_score(fi, ek);
                uint32_t cs = component_size[static_cast<size_t>(face_component[fi])];
                bool better;
                if (coh > best_coh + 1e-9 || coh < best_coh - 1e-9) {
                    better = (coh > best_coh);
                } else if (cs != best_cs) {
                    better = (cs > best_cs);
                } else {
                    better = (fi < best_fi);
                }
                if (better) {
                    best_fi = fi;
                    best_coh = coh;
                    best_cs = cs;
                }
            }
            for (size_t i = 0; i < group.size(); ++i) {
                uint32_t fi = group[i];
                if (fi != best_fi) {
                    remove_face[fi] = true;
                    result.faces_removed++;
                }
            }
        };

        // Less-aggressive policy (2026-05-13, per `algorithm_revision.md`
        // direction #3 and the recent design audit). When BOTH groups have
        // ≥ 2 faces AND the faces within each group form a COHERENT
        // surface (high pairwise normal alignment, avg dot ≥ 0.85), this
        // is a genuine "two surfaces meeting at one edge" pattern (e.g.
        // mug handle/body, trumpet bell/tube, CAD T-junctions). Preserve
        // them — the NM edge stays in topology, but the input geometry
        // is retained.
        //
        // Stripping is still applied:
        //   - when one group has only 1 face (stray fin → remove);
        //   - when intra-group coherence is low (FWN marching-cubes pinch
        //     points: faces with random normals → not legitimate surfaces);
        //   - when fan size is huge (≥ 4 in a single group: pathological).
        auto intra_group_coherence = [&](const std::vector<uint32_t>& g) -> double {
            if (g.size() < 2) return 1.0;
            double sum = 0.0;
            int cnt = 0;
            for (size_t a = 0; a < g.size(); ++a) {
                const internal::Vec3d na = face_normal(g[a]);
                for (size_t b = a + 1; b < g.size(); ++b) {
                    const internal::Vec3d nb = face_normal(g[b]);
                    sum += internal::dot(na, nb);
                    ++cnt;
                }
            }
            return cnt > 0 ? (sum / static_cast<double>(cnt)) : 0.0;
        };
        const bool small_fans   = (groupA.size() <= 3 && groupB.size() <= 3);
        const bool both_multi   = (groupA.size() >= 2 && groupB.size() >= 2);
        const double cohA       = intra_group_coherence(groupA);
        const double cohB       = intra_group_coherence(groupB);
        const bool both_coherent = (cohA >= 0.85 && cohB >= 0.85);
        // Size gate: preserve only when input mesh is LARGE (>= 1000 faces).
        // Real-world CAD/scan inputs with legitimate multi-surface junctions
        // (mug, trumpet, t10k_*) are big; FWN+LevelSet subcomponents fed
        // back through this stage are typically smaller and benefit from
        // the aggressive strip behaviour. Same rationale as the Phase 8R
        // catastrophic-collapse guard.
        const bool large_input  = (nfaces >= 100);
        const bool preserve_multi = both_multi && small_fans && both_coherent && large_input;
        if (!preserve_multi) {
            mark_excess(groupA);
            mark_excess(groupB);
        }
        result.edges_fixed++;
    }

    // Step 5: Build output mesh excluding removed faces, then compact vertices.
    std::vector<bool> used_vert(nverts, false);
    std::vector<std::array<uint32_t, 3>> new_faces;
    new_faces.reserve(nfaces - result.faces_removed);

    for (uint32_t fi = 0; fi < nfaces; ++fi) {
        if (remove_face[fi]) continue;
        const auto& f = mesh.faces[fi];
        new_faces.push_back(f);
        for (int k = 0; k < 3; ++k)
            used_vert[f[k]] = true;
    }

    // Build vertex remap: old index -> new index
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

    // Reindex faces using remap
    for (size_t fi = 0; fi < new_faces.size(); ++fi) {
        auto& f = new_faces[fi];
        for (int k = 0; k < 3; ++k)
            f[k] = remap[f[k]];
    }

    result.mesh.vertices = std::move(new_verts);
    result.mesh.faces = std::move(new_faces);
    return result;
}

} // namespace meshseal::stages
