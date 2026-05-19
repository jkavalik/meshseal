#pragma once

#include "vec3.h"
#include <array>
#include <utility>
#include <vector>

namespace meshseal::internal {

// ────────────────────────────────────────────────────────────────────────
// arrangement.h
//
// Triangle-triangle intersection and triangle-by-cut-segment splitting.
// Used by soup_reconstruct's Step 2 to discover natural meeting lines
// between overlapping/interpenetrating input triangles and split them at
// those lines, so the voxel-oracle classifier can then decide per-fragment
// what survives.
//
// Scope of this v1:
//   - Non-coplanar pairs only. Coplanar pairs return no intersection. (Rare
//     in practice; coplanar handling adds significant code without a
//     proportional fixture-rate win.)
//   - The triangle split assumes cut endpoints lie ON the triangle's edges
//     (with tolerance). Cuts whose endpoints fall in the triangle's interior
//     (e.g. another triangle ending inside this one) are dropped — they
//     would require a full PSLG / constrained-Delaunay implementation that
//     is out of scope here. Most fixtures the soup pipeline cares about
//     have cuts that run edge-to-edge.
// ────────────────────────────────────────────────────────────────────────

// Möller's triangle-triangle intersection.
//
// Computes the line segment along which two non-coplanar triangles
// intersect, clipped to both interiors. Returns:
//   - empty vector if the triangles do not intersect, are coplanar, or are
//     degenerate (zero-area).
//   - two 3D points otherwise (a 3D segment).
std::vector<Vec3d> tri_tri_intersect(
    const Vec3d& a1, const Vec3d& b1, const Vec3d& c1,
    const Vec3d& a2, const Vec3d& b2, const Vec3d& c2);


// Split a triangle by a set of cut segments that lie on (or very near) the
// triangle's plane and have endpoints on the triangle's edges. Cuts are
// applied recursively — each cut splits one existing piece into two,
// pieces are re-fed for further cuts. Order is preserved (cuts[0] applied
// first, then cuts[1], …).
//
// `edge_tol` is the absolute distance below which a cut endpoint is
// considered "on" a triangle edge. Pick relative to your mesh scale (a
// fraction of mean edge length works well).
//
// Returns a list of triangles whose union equals the original triangle
// (up to fp error introduced by cut-edge interpolation). If a cut cannot
// be applied to a piece (endpoint not on any edge of that piece, both
// endpoints on the same edge, etc.) it is silently skipped for that piece.
std::vector<std::array<Vec3d, 3>> split_triangle_by_cuts(
    const Vec3d& v0, const Vec3d& v1, const Vec3d& v2,
    const std::vector<std::pair<Vec3d, Vec3d>>& cuts,
    double edge_tol);

} // namespace meshseal::internal
