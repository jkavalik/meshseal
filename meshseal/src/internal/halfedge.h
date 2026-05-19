#pragma once
#include "../../include/meshseal/meshseal.h"
#include <vector>
#include <unordered_map>
#include <cstdint>

namespace meshseal::internal {

// Half-edge structure built from a Mesh.
// For each directed edge (u→v), stores the incident face index.
// "twin" edges are (v→u).
struct HalfEdge {
    uint32_t vertex_from;
    uint32_t vertex_to;
    uint32_t face;        // index into mesh.faces
    int32_t  twin = -1;  // index of the twin half-edge, -1 if boundary
    uint32_t next;        // next half-edge in the same face (CCW)
    uint32_t prev;        // prev half-edge in the same face
};

struct HalfEdgeMesh {
    std::vector<HalfEdge> half_edges;

    // Map from directed edge (from, to) to half-edge index
    std::unordered_map<uint64_t, uint32_t> edge_map;

    // For each face: the index of its first half-edge (the one from face[0]→face[1])
    std::vector<uint32_t> face_start;

    // Build from mesh. Requires no degenerate or duplicate faces.
    static HalfEdgeMesh build(const Mesh& mesh);

    // Encode a directed edge as a uint64_t key
    static uint64_t edge_key(uint32_t from, uint32_t to) {
        return (static_cast<uint64_t>(from) << 32) | static_cast<uint64_t>(to);
    }

    // Returns true if edge (from, to) exists in the half-edge structure
    bool has_edge(uint32_t from, uint32_t to) const {
        return edge_map.count(edge_key(from, to)) > 0;
    }

    // Returns the half-edge index for directed edge (from, to), or -1 if not found
    int32_t find_half_edge(uint32_t from, uint32_t to) const {
        auto it = edge_map.find(edge_key(from, to));
        return (it != edge_map.end()) ? static_cast<int32_t>(it->second) : -1;
    }
};

} // namespace meshseal::internal
