#pragma once
#include "../../include/meshseal/meshseal.h"
#include <cstdint>

namespace meshseal::stages {

struct CoplanarFanDropResult {
    Mesh     mesh;
    bool     applied      = false;
    uint32_t fans_dropped = 0;   // count of resolved 4-fans
    uint32_t faces_removed = 0;  // total faces removed (always 2 per fan)
};

// Resolve same-sided coplanar 4-fans on non-manifold edges.
//
// Pattern: a non-manifold edge with exactly 4 incident faces, all in the
// same plane (|n·n'| > 0.99 for every pair), with the edge appearing
// twice in each winding direction (2 forward + 2 reverse). This is
// signature of a doubled coplanar surface — two overlapping triangulations
// of the same plane region sharing one chord. Common on multi-color 3MF
// exports (paint-on-part) where adjacent color regions overlap at their
// boundaries, and on CSG-union outputs where coplanar regions met at an
// edge.
//
// Resolution: identify the two manifold pairings (each forward face goes
// with one reverse face); the pairing whose total face-area is smaller
// is the "duplicate layer" — drop those two faces. The remaining 2
// faces give a clean 2-manifold pair on the edge. The drop opens 8 new
// boundary edges (each dropped face had 2 outer edges that connected to
// surrounding manifold geometry); the caller is responsible for running
// fill_holes or recursive repair() afterwards to close them.
//
// This stage is structural — not a heuristic local repair. The
// invariant is: no face is dropped unless it's part of a strict 2+2
// coplanar fan around an NM edge. Single-fan-pair cases, antiparallel
// cases (`coplanar_dedup` target), and >4-incidence cases are all
// untouched.
CoplanarFanDropResult coplanar_fan_drop(const Mesh& mesh,
                                        double coplanar_dot = 0.99);

} // namespace meshseal::stages
