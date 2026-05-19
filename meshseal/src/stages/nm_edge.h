#pragma once
#include "../../include/meshseal/meshseal.h"
#include <cstdint>

namespace meshseal::stages {

struct NmEdgeResult {
    Mesh mesh;
    uint32_t faces_removed;
    uint32_t edges_fixed;          // how many non-manifold edges were resolved
    bool had_non_manifold_edges;   // true if any nm edges were found
};

// Remove faces that cause non-manifold edges (edges shared by >= 3 faces).
// Strategy: for each group of faces sharing an undirected edge, keep only 2 faces
// (the pair with opposite winding, i.e., consistent manifold orientation). If no
// consistent pair exists, keep only the 2 faces that belong to the largest connected
// component. Excess faces are "fins" and are removed.
NmEdgeResult fix_non_manifold_edges(const Mesh& mesh);

} // namespace meshseal::stages
