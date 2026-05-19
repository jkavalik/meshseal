#pragma once
#include "../../include/meshseal/meshseal.h"
#include <cstdint>
#include <vector>

namespace meshseal::internal {

// Classification for a connected face-component after Phase 4 (nm_edge).
enum class ComponentClass {
    SOUP,              // Too small / all-boundary / flat+mostly-open  → reconstruction or discard
    OPEN,              // Has boundary edges with 3D extent           → OPEN queue
    NO_BOUNDARY,       // No boundary edges                           → CLOSED queue directly
    NEEDS_RECONSTRUCTION, // Explicitly tagged for reconstruction (Phase 4b output)
};

struct ComponentInfo {
    std::vector<uint32_t> face_indices; // Indices into the source mesh face array
    ComponentClass cls;
    double open_ratio;      // boundary_edge_count / total_edge_count
    double planarity_ratio; // λ₃/λ₁ from PCA of vertex positions (0 = flat, 1 = isotropic)
    uint32_t face_count;
    bool reconstruction_attempted = false; // Set true once soup_reconstruct has been tried
};

// Decompose `mesh` into connected components via edge-connectivity BFS (any shared
// edge, including boundary edges).  Classify each component according to the routing
// rules from the Revised Pipeline Architecture (Phase 5R).
//
// Routing rules (first match wins):
//   face_count ≤ 3                              → SOUP
//   open_ratio == 1.0                           → SOUP
//   planarity_ratio < planarity_thresh
//     AND open_ratio > open_ratio_thresh        → SOUP
//   open_ratio > 0                              → OPEN
//   open_ratio == 0                             → NO_BOUNDARY
//
// planarity_thresh default: 0.01  (configurable)
// open_ratio_thresh default: 0.50 (configurable)
std::vector<ComponentInfo> classify_components(
    const Mesh& mesh,
    double planarity_thresh = 0.01,
    double open_ratio_thresh = 0.50);

// Extract a single component into its own Mesh (compact vertex array, new face indices).
Mesh extract_component(const Mesh& src, const std::vector<uint32_t>& face_indices);

// Concatenate several meshes into one.  Each input mesh gets its own vertex block;
// face indices are adjusted by the per-mesh vertex offset.
Mesh merge_meshes(const std::vector<Mesh>& meshes);

} // namespace meshseal::internal
