#pragma once
#include "../../include/meshseal/meshseal.h"

namespace meshseal::stages {

struct DegenerateResult {
    Mesh mesh;
    uint32_t degenerate_removed;
    uint32_t duplicate_removed;
    uint32_t isolated_vertices_removed;
};

// Remove degenerate faces (zero-area) and duplicate faces.
// Also compacts the vertex list to remove unreferenced (isolated) vertices.
//
// area_threshold: a face is degenerate if its area is strictly below
// `area_threshold`. The default value `0.0` is a sentinel that triggers
// scale-aware threshold derivation inside `remove_degenerate_faces`
// (currently 1e-22 × bbox_diag²), which catches truly-collinear zero-area
// triangles without filtering legitimate sliver triangles produced by
// hole-fill on real-world meshes. Pass an explicit positive value to use
// an absolute threshold instead.
DegenerateResult remove_degenerate(const Mesh& mesh, double area_threshold = 1e-14);

// Remove only zero-area faces (no duplicate check, no vertex compaction).
// See `remove_degenerate` for the `area_threshold` semantics.
Mesh remove_degenerate_faces(const Mesh& mesh, double area_threshold = 1e-14);

// Remove only duplicate faces (exact duplicate by sorted vertex index tuple).
// Does not compact vertices.
Mesh remove_duplicate_faces(const Mesh& mesh);

} // namespace meshseal::stages
