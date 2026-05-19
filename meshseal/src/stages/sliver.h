#pragma once
#include "../../include/meshseal/meshseal.h"
#include <cstdint>

namespace meshseal::stages {

struct SliverResult {
    Mesh     mesh;
    uint32_t slivers_collapsed = 0;  // number of short-edge collapses applied
};

// Collapse sliver triangles by short-edge collapse.
//
// A "sliver" is a triangle with very low shape quality
//   q = 4·√3·area / (sum of squared edge lengths)
// (Pillinger's measure: 1 for an equilateral triangle, → 0 for a sliver).
// Extreme needles (q below `quality_threshold`) are the target — e.g. the
// thin-strip flap faces left by imperfect CSG self-intersection resolution
// (t10k_1582375: flap faces are q≈0.001, edges [0.004, 6, 6]).
//
// For each sliver, the shortest edge is collapsed (its two endpoints merged)
// IF that edge is manifold (exactly 2 incident faces). Collapsing the short
// edge removes the needle without opening a boundary. Non-manifold or
// boundary shortest edges are skipped — collapsing those is unsafe (it would
// tear the surface), so a sliver whose short edge IS the non-manifold edge
// is left untouched.
//
// Guards (a collapse is rejected if any fails):
//   * normal-flip: no surviving face incident to the merged vertex may flip
//     its normal (dot of pre/post normal must stay positive);
//   * non-manifold creation: the collapse must not push any edge to ≥3 faces;
//   * degenerate creation: the collapse must not create a new zero-area face
//     other than the two faces of the collapsed edge (which are removed).
//
// Geometry preservation: the collapsed edge is short by construction, so the
// merged vertex moves at most one short-edge length — sub-pixel for the
// needle flaps this targets.
SliverResult collapse_slivers(const Mesh& mesh,
                              double quality_threshold = 0.02);

} // namespace meshseal::stages
