#pragma once
#include "../../include/meshseal/meshseal.h"
#include <cstdint>

namespace meshseal::stages {

struct ThinResult {
    Mesh mesh;
    uint32_t faces_removed; // thin faces removed
    uint32_t pairs_found;   // number of back-to-back face pairs detected
};

// Remove back-to-back face pairs (thin features / fins).
// Two faces are considered a thin pair if:
//   1. Their centroids are within `distance_threshold` of each other.
//   2. Their normals have dot product < -cos_threshold (facing roughly opposite).
//      Default cos_threshold = cos(150°) = -0.866 (normals within 30° of antiparallel).
// Both faces in each pair are removed.
ThinResult remove_thin_features(const Mesh& mesh,
                                double distance_threshold = 1e-3,
                                double cos_threshold = -0.866);

} // namespace meshseal::stages
