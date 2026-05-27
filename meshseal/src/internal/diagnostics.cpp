#include "diagnostics.h"
#include "vec3.h"
#include "../../include/meshseal/meshseal.h"

#include <algorithm>
#include <array>
#include <map>
#include <queue>
#include <set>
#include <unordered_map>
#include <vector>

namespace meshseal::internal {

MeshDiagnostics compute_diagnostics(const meshseal::Mesh& mesh) {
    MeshDiagnostics diag;

    const auto& verts = mesh.vertices;
    const auto& faces = mesh.faces;

    if (faces.empty()) {
        return diag;
    }

    // --- Edge manifoldness: count how many faces reference each undirected edge ---
    // key: (min_v, max_v)
    std::map<std::pair<uint32_t, uint32_t>, int> edge_count;

    for (const auto& f : faces) {
        for (int e = 0; e < 3; ++e) {
            const uint32_t a = f[e];
            const uint32_t b = f[(e + 1) % 3];
            const auto key = std::make_pair(std::min(a, b), std::max(a, b));
            edge_count[key]++;
        }
    }

    for (const auto& [edge, count] : edge_count) {
        if (count == 1) {
            ++diag.open_boundary_edges;
        } else if (count >= 3) {
            ++diag.non_manifold_edges;
        }
    }

    // is_orientable: conservative — true iff no non-manifold edges
    diag.is_orientable = (diag.non_manifold_edges == 0);

    // --- Degenerate faces ---
    for (const auto& f : faces) {
        const Vec3d a = verts[f[0]];
        const Vec3d b = verts[f[1]];
        const Vec3d c = verts[f[2]];
        if (triangle_area(a, b, c) < 1e-14) {
            ++diag.degenerate_faces;
        }
    }

    // --- Duplicate faces ---
    // Normalise each face to a sorted triple for comparison
    std::map<std::array<uint32_t, 3>, int> face_count;
    for (const auto& f : faces) {
        std::array<uint32_t, 3> sorted = f;
        std::sort(sorted.begin(), sorted.end());
        face_count[sorted]++;
    }
    for (const auto& [key, count] : face_count) {
        if (count > 1) {
            diag.duplicate_faces += count - 1;
        }
    }

    // --- Signed volume ---
    // (signed_tet_volume implicitly uses the origin as the apex)
    for (const auto& f : faces) {
        const Vec3d a = verts[f[0]];
        const Vec3d b = verts[f[1]];
        const Vec3d c = verts[f[2]];
        diag.signed_volume += signed_tet_volume(a, b, c);
    }

    // --- Connected components via BFS on vertex adjacency ---
    const uint32_t n_verts = static_cast<uint32_t>(verts.size());
    std::vector<std::vector<uint32_t>> adj(n_verts);
    for (const auto& f : faces) {
        for (int e = 0; e < 3; ++e) {
            const uint32_t a = f[e];
            const uint32_t b = f[(e + 1) % 3];
            adj[a].push_back(b);
            adj[b].push_back(a);
        }
    }

    // Only visit vertices that are referenced by at least one face
    std::vector<bool> has_face(n_verts, false);
    for (const auto& f : faces) {
        has_face[f[0]] = true;
        has_face[f[1]] = true;
        has_face[f[2]] = true;
    }

    std::vector<bool> visited(n_verts, false);
    int components = 0;
    for (uint32_t start = 0; start < n_verts; ++start) {
        if (!has_face[start] || visited[start]) {
            continue;
        }
        ++components;
        std::queue<uint32_t> q;
        q.push(start);
        visited[start] = true;
        while (!q.empty()) {
            const uint32_t v = q.front();
            q.pop();
            for (const uint32_t nb : adj[v]) {
                if (!visited[nb]) {
                    visited[nb] = true;
                    q.push(nb);
                }
            }
        }
    }

    diag.component_count = components;
    diag.is_connected = (components == 1);

    // non_manifold_vertices and isolated_vertices: placeholder for Phase 4
    diag.non_manifold_vertices = 0;
    diag.isolated_vertices = 0;

    return diag;
}

} // namespace meshseal::internal
