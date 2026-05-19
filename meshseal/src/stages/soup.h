#pragma once
#include "../../include/meshseal/meshseal.h"
#include <cstdint>

namespace meshseal::stages {

struct SoupResult {
    Mesh mesh;
    bool success;          // true if manifold successfully reconstructed the mesh
    bool was_needed;       // true if input was non-manifold (soup reconstruction actually ran)
    uint32_t faces_before;
    uint32_t faces_after;
};

// Last-resort triangle soup reconstruction using the manifold library.
// Input should be a mesh that failed all prior stages.
// Attempts to create a manifold solid via manifold's triangle soup import.
SoupResult reconstruct_soup(const Mesh& mesh);

} // namespace meshseal::stages
