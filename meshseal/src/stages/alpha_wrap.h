#pragma once
#include "../../include/meshseal/meshseal.h"
#include <cstdint>
#include <string>

namespace meshseal::stages {

struct AlphaWrapResult {
    Mesh        mesh;
    bool        success    = false;
    double      alpha_used = 0.0;
    std::string reason;
};

// Alpha-wrapping: reconstruct a watertight surface by morphologically
// CLOSING the input solid with a ball of radius alpha.
//
// A finite-radius probe cannot enter gaps narrower than 2*alpha, so the
// closing bridges the open seams of a structured-soup input (LDraw/LEGO
// exports) instead of leaking through them. Concavities and gaps narrower
// than 2*alpha are sealed; locally-convex boundary is left unchanged.
//
// Algorithm (closed form on a voxel grid — not a per-angle sweep):
//   1. voxelize the bbox; seed Surface voxels by point-sampling triangles;
//   2. distance transform #1 -> dist_to_surface;
//   3. flood-fill from the grid boundary through {dist_to_surface > alpha}
//      -> the reachable probe-centre set C;
//   4. distance transform #2 -> dist_to_C;
//   5. signed field  s = alpha - dist_to_C   (s < 0 inside the solid);
//   6. manifold::LevelSet(s) -> watertight mesh, welded by exact float32.
//
// This is a remesh — it does NOT preserve the original triangles — so it is
// intended only as a gated last-resort reconstructor for inputs the
// geometry-preserving pipeline cannot repair (see repair()'s destruction
// fallback). It cannot solidify an open shell that bounds no volume, cannot
// help a non-orientable input, and fills enclosed cavities.
//
// `alpha` < 0 -> auto (a small fraction of bbox diagonal, floored so the
// ball spans several voxels). `voxel_res` is the grid resolution along the
// longest bbox axis.
AlphaWrapResult alpha_wrap(const Mesh& mesh,
                           double alpha     = -1.0,
                           int    voxel_res = 160);

} // namespace meshseal::stages
