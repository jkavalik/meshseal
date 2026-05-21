#pragma once
#include "../../include/meshseal/meshseal.h"
#include <cstdint>

namespace meshseal::internal {

// Count distinct self-intersecting triangle pairs (Möller T-T).
//
// Pruning: bucket triangle bounding boxes into a uniform grid sized by the
// median triangle edge length; test only pairs that share at least one
// voxel. Shared-edge pairs (adjacent triangles in a manifold mesh) are
// excluded — they intersect along the shared edge by definition.
//
// Returns the number of unordered triangle pairs (i,j) with i<j that
// genuinely cross. Coplanar pairs and shared-edge adjacency don't count.
// O(F * k) where k is the bucket population — effectively O(F) on
// well-distributed meshes.
int count_self_intersections(const Mesh& mesh);

} // namespace meshseal::internal
