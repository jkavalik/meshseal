#pragma once
#include "../../include/meshseal/meshseal.h"
#include <cstdint>

namespace meshseal::stages {

struct BridgeLoopsResult {
    Mesh     mesh;
    uint32_t bridges_made = 0;   // number of loop pairs zipped
    uint32_t faces_added  = 0;   // bridge-strip triangles appended
};

// Bridge pairs of open boundary loops that belong to DIFFERENT connected
// components when their centroids are close relative to their size.
//
// Motivation: a model can have a topological opening shared by two shells —
// e.g. a head shell with a mouth hole and a separate mouth-cavity shell with
// its own opening (captain_toad). The two rims are distinct boundary loops in
// distinct components. Hole-fill, run per component, would cap each loop
// independently — sealing the shared cavity into a flat membrane and
// destroying the open-mouth geometry. The correct repair zips the two rims
// together with a triangle strip, leaving the cavity open.
//
// Algorithm (ported from the meshseal_go predecessor's bridge_paired_loops):
//   1. find edge-connected components;
//   2. trace open boundary loops, keep those with [min,max]_loop_verts;
//   3. form cross-component loop pairs, sort by centroid distance;
//   4. greedily accept a pair when centroid_dist < max_dist_factor * avg_radius
//      and the centroid-to-centroid segment is not blocked by solid geometry
//      (Möller–Trumbore line-of-sight test); each loop bridges at most once;
//   5. zip each accepted pair with a greedy shortest-edge triangle strip.
//
// Runs before component classification / hole-fill. No-op when the mesh has
// <2 components or no qualifying cross-component loop pair.
// `max_dist_factor`, `min_radius_ratio`, `min_vert_ratio` form a tight
// "genuinely-coincident loop pair" gate: a real shared opening has two rims
// that nearly overlap and are similarly sized. Loose pairings (two unrelated
// holes that merely happen to be cross-component and roughly near each other)
// are rejected — they would otherwise zip unrelated geometry together.
BridgeLoopsResult bridge_paired_loops(const Mesh& mesh,
                                      double   max_dist_factor  = 0.35,
                                      double   min_radius_ratio = 0.5,
                                      double   min_vert_ratio   = 0.5,
                                      uint32_t min_loop_verts   = 6,
                                      uint32_t max_loop_verts   = 400);

} // namespace meshseal::stages
