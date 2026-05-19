#pragma once
#include "../../include/meshseal/meshseal.h"
#include <cstdint>

namespace meshseal::stages {

struct IntersectResult {
    Mesh mesh;
    bool had_intersections; // true if manifold detected and fixed intersections
    bool manifold_failed;   // true if manifold library returned error/empty
    uint32_t faces_before;
    uint32_t faces_after;
};

// Resolve self-intersections using manifold boolean self-union.
// Returns the input mesh unchanged if manifold fails or detects no issues.
IntersectResult resolve_intersections(const Mesh& mesh);

} // namespace meshseal::stages
