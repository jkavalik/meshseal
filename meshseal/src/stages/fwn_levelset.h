#pragma once
#include "../../include/meshseal/meshseal.h"
#include <cstdint>
#include <string>

namespace meshseal::stages {

// Fast Winding Number (Jacobson 2013) + manifold::LevelSet marching cubes.
//
// Stage 10. Designed as a last-resort for triangle soups that the existing
// arrangement-based soup_reconstruct cannot close (random orientation soups
// like soup_seed*, isolated triangle clouds, back-to-back face pairs). The
// pipeline calls this when soup_reconstruct returns empty or volume-zero
// output, so the regression-free guarantee is: if FWN+LevelSet produces no
// valid closed mesh, the caller retains the original (failing) output.
//
// Algorithm:
//   1. Build a per-point winding number function ω(p) from the input
//      triangles (Jacobson 2013 closed-form formula, O(F) per query).
//   2. Define SDF s(p) = 0.5 - ω(p) so:
//          inside  ω > 0.5  ⇒  s < 0
//          outside ω < 0.5  ⇒  s > 0
//   3. Call manifold::LevelSet(s, bounds, edge_length).
//
// `samples_per_axis` controls marching-cubes resolution: edge_length =
// bbox_diag / samples_per_axis. 64 is a good default balance.
//
// `wn_threshold` selects the isosurface (default 0.5). Lower values
// (e.g. 0.3) include more geometry, useful for sparse soups; higher
// values (e.g. 0.7) are stricter.
struct FwnLevelSetResult {
    Mesh        mesh;
    bool        success = false;
    std::string reason;
    uint32_t    faces_in  = 0;
    uint32_t    faces_out = 0;
};

FwnLevelSetResult fwn_levelset(const Mesh& mesh,
                               int    samples_per_axis = 20,
                               double wn_threshold     = 0.5);

// Voxel-occupancy + manifold::LevelSet — the soup reconstruction the original
// design docs prescribed and Stage 10 should have been doing all along.
//
// Algorithm (per repair_algorithms_investigation.md §1.3, the random-soup
// branch the originals mandated):
//   1. Rasterize all input triangles into a VoxelGrid → Surface cells.
//   2. 6-connected flood fill from grid boundary → Outside cells.
//   3. Cells still unmarked → Inside (an enclosed pocket of the input).
//   4. Define SDF s(p) such that s<0 inside, s>0 outside, s=0 on surface,
//      sampled from the voxel grid. Feed to manifold::LevelSet.
//
// Orientation-agnostic (unlike FWN, which requires consistent winding to
// produce a meaningful winding number). This is what unlocks the
// soup_seed* family where FWN's winding number is garbage on random
// triangle orientations.
//
// `voxel_res` is grid resolution along the longest bbox axis (default 64);
// `samples_per_axis` controls marching-cubes density (default 20). They can
// differ — the voxel oracle wants finer resolution than the output mesh.
FwnLevelSetResult voxel_levelset(const Mesh& mesh,
                                 int voxel_res        = 64,
                                 int samples_per_axis = 20);

} // namespace meshseal::stages
