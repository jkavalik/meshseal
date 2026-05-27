#include "halfedge.h"
#include <algorithm>

namespace meshseal::internal {

HalfEdgeMesh HalfEdgeMesh::build(const Mesh& mesh) {
    HalfEdgeMesh hm;
    const uint32_t nfaces = static_cast<uint32_t>(mesh.faces.size());

    hm.half_edges.resize(nfaces * 3);
    hm.face_start.resize(nfaces);

    // Step 1: create half-edges and populate edge_map.
    // Pre-condition: input mesh has no duplicate directed edges
    // (no two faces traverse the same edge in the same direction).
    // This is implicit in "manifold + correctly oriented" but easy
    // to violate on raw triangle soup. The map insert below uses
    // `try_emplace`; collisions track in `hm.directed_edge_collisions`
    // so a future caller can detect and reject the input rather than
    // silently use the first-or-last write.
    hm.directed_edge_collisions = 0;
    for (uint32_t i = 0; i < nfaces; ++i) {
        const auto& f = mesh.faces[i];
        const uint32_t v0 = f[0], v1 = f[1], v2 = f[2];

        const uint32_t he0 = 3 * i + 0;
        const uint32_t he1 = 3 * i + 1;
        const uint32_t he2 = 3 * i + 2;

        hm.half_edges[he0] = {v0, v1, i, -1, he1, he2};
        hm.half_edges[he1] = {v1, v2, i, -1, he2, he0};
        hm.half_edges[he2] = {v2, v0, i, -1, he0, he1};

        const std::pair<uint64_t, uint32_t> kvs[3] = {
            {edge_key(v0, v1), he0},
            {edge_key(v1, v2), he1},
            {edge_key(v2, v0), he2},
        };
        for (const auto& kv : kvs) {
            auto inserted = hm.edge_map.try_emplace(kv.first, kv.second);
            if (!inserted.second) ++hm.directed_edge_collisions;
        }

        hm.face_start[i] = he0;
    }

    // Step 2: link twins
    for (uint32_t idx = 0; idx < static_cast<uint32_t>(hm.half_edges.size()); ++idx) {
        HalfEdge& he = hm.half_edges[idx];
        auto it = hm.edge_map.find(edge_key(he.vertex_to, he.vertex_from));
        if (it != hm.edge_map.end()) {
            he.twin = static_cast<int32_t>(it->second);
        }
    }

    return hm;
}

} // namespace meshseal::internal
