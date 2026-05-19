#pragma once
#include "../../include/meshseal/meshseal.h"
#include <cstdint>

namespace meshseal::stages {

struct HoleFillResult {
    Mesh mesh;
    uint32_t holes_filled;
    uint32_t faces_added;
};

// Fill boundary loops. Small loops (≤ max_fan_size): fan from centroid. Larger: ear-clipping.
HoleFillResult fill_holes(const Mesh& mesh, uint32_t max_fan_size = 20);

} // namespace meshseal::stages
