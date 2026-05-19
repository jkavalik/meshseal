#pragma once
#include "../../include/meshseal/meshseal.h"
#include <cstdint>
#include <vector>

namespace meshseal::stages {

struct ShellInfo {
    uint32_t component_index; // 0-based
    uint32_t face_count;
    double   volume;          // signed; positive = outward normals
    bool     dropped;         // true if removed (below volume threshold)
    bool     is_contained;    // true if this shell is inside another shell
};

struct ShellResult {
    Mesh mesh;
    std::vector<ShellInfo> shells;
    uint32_t shells_dropped;
    uint32_t shells_kept;
};

// Analyze and filter shells.
// vol_threshold: drop shells with |vol| < vol_threshold * max(|vol|).
// Default lowered to 0.001 (was 0.01) — the 1 % threshold was dropping small
// but legitimate features (handles, buttons, tiny components of complex
// real-world inputs) and contributing to Hausdorff regressions on
// preservation-oriented inputs.
ShellResult analyze_shells(const Mesh& mesh, double vol_threshold = 0.001);

} // namespace meshseal::stages
