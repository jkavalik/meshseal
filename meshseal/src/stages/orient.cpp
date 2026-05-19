#include "orient.h"
#include "../internal/vec3.h"
#include <algorithm>
#include <cmath>
#include <queue>
#include <unordered_map>
#include <vector>

namespace meshseal::stages {

namespace {

// Pack an undirected edge as key: min first
static uint64_t undirected_edge_key(uint32_t a, uint32_t b) {
    if (a > b) std::swap(a, b);
    return (static_cast<uint64_t>(a) << 32) | static_cast<uint64_t>(b);
}

} // anonymous namespace

OrientResult orient_mesh(const Mesh& mesh, bool do_signed_volume_flip) {
    OrientResult res;
    res.faces_flipped = 0;
    res.components_flipped = 0;
    res.was_orientable = true;

    const uint32_t nfaces = static_cast<uint32_t>(mesh.faces.size());
    if (nfaces == 0) {
        res.mesh = mesh;
        return res;
    }

    // Working copy of faces (winding may be adjusted)
    std::vector<std::array<uint32_t, 3>> faces(mesh.faces.begin(), mesh.faces.end());

    // Step 1: build undirected edge → faces map
    // undirected_edge_to_faces: edge → list of face indices sharing it
    std::unordered_map<uint64_t, std::vector<uint32_t>> undirected_edge_to_faces;
    undirected_edge_to_faces.reserve(nfaces * 3);

    for (uint32_t i = 0; i < nfaces; ++i) {
        const auto& f = faces[i];
        const uint32_t verts[3] = {f[0], f[1], f[2]};
        for (int e = 0; e < 3; ++e) {
            const uint32_t u = verts[e];
            const uint32_t v = verts[(e + 1) % 3];
            undirected_edge_to_faces[undirected_edge_key(u, v)].push_back(i);
        }
    }

    // Detect non-manifold edges (shared by 3+ faces) and mark them
    // We will skip BFS propagation across such edges
    std::unordered_map<uint64_t, bool> non_manifold_edge;
    for (auto it = undirected_edge_to_faces.begin(); it != undirected_edge_to_faces.end(); ++it) {
        if (it->second.size() > 2) {
            non_manifold_edge[it->first] = true;
            res.was_orientable = false;
        }
    }

    // Build face adjacency for BFS: face_neighbors[i] = list of (neighbor_face, shared_edge_u, shared_edge_v)
    // where (u,v) is the directed edge as it appears in face i
    struct FaceEdge {
        uint32_t neighbor;
        uint32_t u, v; // directed edge in face i
    };
    std::vector<std::vector<FaceEdge>> face_neighbors(nfaces);
    for (uint32_t i = 0; i < nfaces; ++i) {
        const auto& f = faces[i];
        const uint32_t verts[3] = {f[0], f[1], f[2]};
        for (int e = 0; e < 3; ++e) {
            const uint32_t u = verts[e];
            const uint32_t v = verts[(e + 1) % 3];
            const uint64_t ukey = undirected_edge_key(u, v);
            // skip non-manifold edges
            if (non_manifold_edge.count(ukey)) continue;
            const auto& shared = undirected_edge_to_faces.at(ukey);
            for (uint32_t nb : shared) {
                if (nb != i) {
                    face_neighbors[i].push_back({nb, u, v});
                }
            }
        }
    }

    // Step 2: BFS over connected components
    std::vector<int> component(nfaces, -1);
    int num_components = 0;
    std::vector<bool> flipped(nfaces, false);

    for (uint32_t seed = 0; seed < nfaces; ++seed) {
        if (component[seed] != -1) continue;

        const int comp_id = num_components++;
        std::queue<uint32_t> q;
        q.push(seed);
        component[seed] = comp_id;

        while (!q.empty()) {
            const uint32_t cur = q.front();
            q.pop();

            for (const auto& fe : face_neighbors[cur]) {
                const uint32_t nb = fe.neighbor;
                if (component[nb] != -1) continue;

                component[nb] = comp_id;

                // Determine actual edge direction in cur by inspecting the current face winding.
                // (Using directed_edge_to_face would be unreliable: when two inconsistently-wound
                // faces share the same directed edge u→v, the map only stores the last writer,
                // causing false "consistent" judgments.)
                const auto& fcur = faces[cur];
                bool cur_has_uv = false;
                for (int e = 0; e < 3; ++e) {
                    if (fcur[e] == fe.u && fcur[(e + 1) % 3] == fe.v) {
                        cur_has_uv = true;
                        break;
                    }
                }
                const uint32_t real_u = cur_has_uv ? fe.u : fe.v;
                const uint32_t real_v = cur_has_uv ? fe.v : fe.u;

                // For consistent winding nb should traverse the same edge as real_v → real_u.
                const auto& fnb = faces[nb];
                bool twin_has_vu = false;
                for (int e = 0; e < 3; ++e) {
                    if (fnb[e] == real_v && fnb[(e + 1) % 3] == real_u) {
                        twin_has_vu = true;
                        break;
                    }
                }

                if (!twin_has_vu) {
                    std::swap(faces[nb][1], faces[nb][2]);
                    flipped[nb] = true;
                    res.faces_flipped++;
                }
                q.push(nb);
            }
        }
    }

    // Step 3: compute signed volume per component and flip if negative.
    // Skipped when `do_signed_volume_flip` is false — caller is re-running
    // after a topology change and wants only the BFS winding propagation,
    // not a second per-component sign decision (which can wrongly invert
    // small components whose signed-volume sign changed for unrelated
    // reasons; see repair.cpp re-orient-after-nm_edge revert notes).
    if (do_signed_volume_flip) {
        std::vector<double> comp_volume(num_components, 0.0);
        std::vector<std::vector<uint32_t>> comp_faces(num_components);
        for (uint32_t i = 0; i < nfaces; ++i) {
            comp_faces[component[i]].push_back(i);
        }

        for (int c = 0; c < num_components; ++c) {
            for (uint32_t fi : comp_faces[c]) {
                const auto& f = faces[fi];
                const auto& v0 = mesh.vertices[f[0]];
                const auto& v1 = mesh.vertices[f[1]];
                const auto& v2 = mesh.vertices[f[2]];
                comp_volume[c] += meshseal::internal::signed_tet_volume(
                    meshseal::internal::Vec3d{v0[0], v0[1], v0[2]},
                    meshseal::internal::Vec3d{v1[0], v1[1], v1[2]},
                    meshseal::internal::Vec3d{v2[0], v2[1], v2[2]});
            }
            if (comp_volume[c] < 0.0) {
                res.components_flipped++;
                for (uint32_t fi : comp_faces[c]) {
                    std::swap(faces[fi][1], faces[fi][2]);
                    if (flipped[fi]) {
                        // double-flip cancels out
                        flipped[fi] = false;
                        res.faces_flipped--;
                    } else {
                        flipped[fi] = true;
                        res.faces_flipped++;
                    }
                }
            }
        }
    }

    // Build output mesh
    res.mesh.vertices = mesh.vertices;
    res.mesh.faces.assign(faces.begin(), faces.end());

    return res;
}

} // namespace meshseal::stages
