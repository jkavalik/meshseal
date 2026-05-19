#pragma once
#include "../../include/meshseal/meshseal.h"
#include <cstdint>

namespace meshseal::stages {

struct NmVertexResult {
    Mesh mesh;
    uint32_t vertices_added;  // how many vertex copies were inserted
    uint32_t vertices_split;  // how many original vertices were identified as non-manifold
};

// Split non-manifold vertices so that each vertex's incident faces form a single fan (disk).
// Non-manifold vertices are detected by BFS/DFS over the one-ring: faces are grouped
// into maximal connected fans (connected = share an edge with another face in the ring).
// If a vertex has more than one fan, extra fans get copies of the vertex.
NmVertexResult split_non_manifold_vertices(const Mesh& mesh);

} // namespace meshseal::stages
