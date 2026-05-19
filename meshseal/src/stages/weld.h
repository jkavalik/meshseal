#pragma once
#include "../../include/meshseal/meshseal.h"
#include <optional>

namespace meshseal::stages {

struct WeldResult {
    Mesh mesh;
    uint32_t vertices_before;
    uint32_t vertices_after;
    double tolerance_used;
};

// Weld vertices within tolerance.
// If tolerance is nullopt, auto-compute from mesh scale:
//   tolerance = max(1e-9, mean_edge_length * 1e-4, bbox_diagonal * 1e-7)
// Returns the welded mesh plus statistics.
// The face list is compacted to remove any faces that degenerate to a
// point or line after welding (all 3 indices the same, or any 2 the same).
WeldResult weld(const Mesh& mesh, std::optional<double> tolerance = std::nullopt);

} // namespace meshseal::stages
