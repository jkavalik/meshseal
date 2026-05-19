#include "nm_vertex.h"
#include <algorithm>
#include <cmath>
#include <queue>
#include <unordered_map>
#include <vector>

namespace meshseal::stages {

NmVertexResult split_non_manifold_vertices(const Mesh& mesh) {
    NmVertexResult result;
    result.mesh = mesh;
    result.vertices_added = 0;
    result.vertices_split = 0;

    const uint32_t n_verts = static_cast<uint32_t>(result.mesh.vertices.size());
    const uint32_t n_faces = static_cast<uint32_t>(result.mesh.faces.size());

    if (n_verts == 0 || n_faces == 0)
        return result;

    // Build vertex → incident face indices
    std::vector<std::vector<uint32_t>> vertex_faces(n_verts);
    for (uint32_t fi = 0; fi < n_faces; ++fi) {
        const auto& f = result.mesh.faces[fi];
        for (int k = 0; k < 3; ++k)
            vertex_faces[f[k]].push_back(fi);
    }

    // Process each original vertex only; new vertices >= n_verts are manifold by construction
    for (uint32_t v = 0; v < n_verts; ++v) {
        // Copy the ring — we may call vertex_faces.emplace_back() below which
        // would reallocate the outer vector and invalidate any reference into it.
        const std::vector<uint32_t> ring = vertex_faces[v];
        const uint32_t ring_size = static_cast<uint32_t>(ring.size());
        if (ring_size <= 1)
            continue;

        // Map global face index → local index within the ring
        std::unordered_map<uint32_t, uint32_t> face_to_local;
        face_to_local.reserve(ring_size);
        for (uint32_t i = 0; i < ring_size; ++i)
            face_to_local[ring[i]] = i;

        // Build local adjacency: two ring-faces are adjacent iff they share a vertex != v
        std::vector<std::vector<uint32_t>> adj(ring_size);
        for (uint32_t i = 0; i < ring_size; ++i) {
            const auto& fa = result.mesh.faces[ring[i]];
            for (int ka = 0; ka < 3; ++ka) {
                uint32_t ov = fa[ka];
                if (ov == v) continue;
                // Every ring-face that also contains ov is adjacent to ring-face i
                for (uint32_t fj : vertex_faces[ov]) {
                    if (fj == ring[i]) continue;
                    auto it = face_to_local.find(fj);
                    if (it != face_to_local.end()) {
                        adj[i].push_back(it->second);
                    }
                }
            }
        }

        // BFS: find connected components of the local adjacency graph
        std::vector<int> component(ring_size, -1);
        int num_components = 0;
        for (uint32_t start = 0; start < ring_size; ++start) {
            if (component[start] != -1) continue;
            std::queue<uint32_t> q;
            q.push(start);
            component[start] = num_components;
            while (!q.empty()) {
                uint32_t cur = q.front();
                q.pop();
                for (uint32_t nb : adj[cur]) {
                    if (component[nb] == -1) {
                        component[nb] = num_components;
                        q.push(nb);
                    }
                }
            }
            ++num_components;
        }

        if (num_components <= 1)
            continue;

        // Non-manifold vertex: first component keeps v in place, each extra
        // component gets a NEW vertex offset slightly toward that fan's
        // centroid. The offset is essential — copying v's position verbatim
        // produces two vertices at bit-identical positions; STL output is
        // float32, so on read-back they collapse and the originally-split
        // fans converge into a high-valence non-manifold edge again. The
        // offset is large enough to survive float32 round-trip but small
        // enough to be visually imperceptible (~1 ppm of the local fan
        // radius), and it points along a geometrically-meaningful direction
        // (toward each fan's interior).
        ++result.vertices_split;
        for (int c = 1; c < num_components; ++c) {
            // Centroid of the OTHER vertices of this fan's faces.
            double cen[3] = {0.0, 0.0, 0.0};
            int cnt = 0;
            for (uint32_t i = 0; i < ring_size; ++i) {
                if (component[i] != c) continue;
                const auto& f = result.mesh.faces[ring[i]];
                for (int k = 0; k < 3; ++k) {
                    if (f[k] == v) continue;
                    const auto& fv = result.mesh.vertices[f[k]];
                    cen[0] += fv[0]; cen[1] += fv[1]; cen[2] += fv[2];
                    ++cnt;
                }
            }
            const auto& vp = result.mesh.vertices[v];
            std::array<double, 3> new_pos = vp;
            if (cnt > 0) {
                cen[0] /= cnt; cen[1] /= cnt; cen[2] /= cnt;
                double dir[3] = { cen[0] - vp[0], cen[1] - vp[1], cen[2] - vp[2] };
                double mag = std::sqrt(dir[0]*dir[0] + dir[1]*dir[1] + dir[2]*dir[2]);
                if (mag > 0.0) {
                    // Offset = local fan radius * 1e-6. Survives float32
                    // round-trip (ULP near unit values is ~1e-7) without
                    // visibly perturbing the geometry.
                    const double off = mag * 1e-6;
                    const double inv = 1.0 / mag;
                    new_pos[0] = vp[0] + dir[0] * inv * off;
                    new_pos[1] = vp[1] + dir[1] * inv * off;
                    new_pos[2] = vp[2] + dir[2] * inv * off;
                }
            }
            const uint32_t new_v = static_cast<uint32_t>(result.mesh.vertices.size());
            result.mesh.vertices.push_back(new_pos);
            // Keep vertex_faces in sync so subsequent vertices can look up new_v safely
            vertex_faces.emplace_back();
            ++result.vertices_added;

            for (uint32_t i = 0; i < ring_size; ++i) {
                if (component[i] != c) continue;
                auto& f = result.mesh.faces[ring[i]];
                for (int k = 0; k < 3; ++k) {
                    if (f[k] == v) f[k] = new_v;
                }
                vertex_faces[new_v].push_back(ring[i]);
            }
        }
    }

    return result;
}

} // namespace meshseal::stages
