#pragma once
#include "../../include/meshseal/meshseal.h"
#include <cstdint>

namespace meshseal::stages {

struct SiSplitResult {
    Mesh     mesh;
    uint32_t pairs_found  = 0;   // self-intersecting triangle pairs detected
    uint32_t faces_before = 0;
    uint32_t faces_after  = 0;   // = pairs_found + faces_before (approx)
};

// Resolve self-intersections by re-triangulating crossing triangles.
//
// For each pair of triangles that cross (Möller T-T), the intersection
// segment is recorded on BOTH triangles. After collecting all cuts per
// triangle, co-linear cuts are merged (the same physical chord may be
// detected by Möller multiple times, once per pair), and the triangle is
// split by `split_triangle_by_cuts` into sub-triangles that share the cut
// edges. After this, no two output triangles physically cross — they may
// coincide or meet edge-to-edge, but they do not penetrate.
//
// This is the standard "arrangement-based" approach used by CGAL's
// `corefine_and_compute_union`. meshseal already uses it inside
// `soup_reconstruct` as a preprocessor; this stage exposes it on its own
// so the SI-cleanup is independent of soup reconstruction.
//
// Faces sharing one or more vertices are excluded from the pair test (their
// "intersection" is at the shared vertex/edge and is a numerical artifact,
// not a genuine SI — same convention as `count_self_intersections`).
//
// Caller is responsible for any post-split cleanup (the back-to-back
// duplicate pairs that emerge from overlapping coincident areas, the
// slivers that may be produced near edge-grazing cuts, etc.). The existing
// `thin_features` / `nm_local_repair` / `collapse_slivers` stages handle
// these.
SiSplitResult resolve_self_intersections(const Mesh& mesh);

} // namespace meshseal::stages
