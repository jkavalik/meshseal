#pragma once
#include "../../include/meshseal/meshseal.h"
#include <array>
#include <cstdint>

namespace meshseal::internal {

// Signed solid angle subtended by triangle (a, b, c) viewed from point p,
// in steradians. Sign matches the triangle winding (right-hand rule).
// Jacobson, Kavan, Sorkine-Hornung (SIGGRAPH 2013) robust atan2 form.
//
// Returns 0 when p coincides with a vertex (degenerate).
double tri_solid_angle(const std::array<double, 3>& a,
                       const std::array<double, 3>& b,
                       const std::array<double, 3>& c,
                       const std::array<double, 3>& p);

// Generalized winding number of point p w.r.t. the mesh. For a closed
// orientable surface, ω(p) is exactly +1 inside, 0 outside. For
// open/soup/multi-shell input, ω is fractional but the level set at 0.5
// is a reasonable "best-guess inside" boundary.
//
// Cost: O(F) per query. For per-face orientation correction this is
// O(F²) overall; acceptable up to ~50k F. For larger inputs a BVH-
// accelerated implementation would be needed (deferred).
double generalized_winding_number(const Mesh& mesh,
                                  const std::array<double, 3>& p);

} // namespace meshseal::internal
