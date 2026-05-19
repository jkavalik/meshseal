#pragma once

#include "../../include/meshseal/meshseal.h"

namespace meshseal::internal {

struct MeshDiagnostics {
    int open_boundary_edges = 0;   // edges with only 1 adjacent face
    int non_manifold_edges = 0;    // edges with 3+ adjacent faces
    int non_manifold_vertices = 0; // vertices with disconnected fan (placeholder)
    int degenerate_faces = 0;      // zero-area triangles
    int duplicate_faces = 0;       // duplicate face tuples
    int isolated_vertices = 0;     // vertices with degree 0 (placeholder)
    bool is_orientable = false;
    bool is_connected = false;
    int component_count = 0;
    double signed_volume = 0.0;    // sum of signed tet volumes
};

// Compute diagnostics for a mesh. O(F log F).
MeshDiagnostics compute_diagnostics(const meshseal::Mesh& mesh);

} // namespace meshseal::internal
