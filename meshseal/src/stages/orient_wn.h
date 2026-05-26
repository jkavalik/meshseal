#pragma once
#include "../../include/meshseal/meshseal.h"
#include <cstdint>

namespace meshseal::stages {

struct OrientWnResult {
    Mesh     mesh;
    bool     applied         = false;
    uint32_t faces_flipped   = 0;
    uint32_t antipar_before  = 0;
    uint32_t antipar_after   = 0;
};

// Per-face orientation correction using the generalized winding number.
//
// Target pattern: INVERSION-type antipar manifold-edge pairs - a
// manifold edge where the two incident faces traverse the edge IN THE
// SAME direction (inconsistent winding) and consequently have anti-
// parallel normals. This is a true orientation defect (a face was
// authored with reversed winding); flipping the inverted faces
// restores consistent outward orientation.
//
// NOT the target: FOLD-type antipar manifold-edge pairs - a manifold
// edge where the two incident faces traverse it in OPPOSITE directions
// (consistent winding) but with anti-parallel normals (because the
// surface folds back on itself at the edge, like a thin pressed-together
// flap). These are membranes / internal flaps, not orientation defects.
// They require coplanar_dedup territory (2D CSG, mismatched-
// triangulation merging). The Bee_v3.3mf vol[3] case is THIS pattern,
// not the inversion pattern - orient_wn correctly leaves it untouched.
//
// For each face F, sample two points: p_out = centroid + eps*N and
// p_in = centroid - eps*N. Compute GWN at each. For a correctly-
// oriented face, GWN(p_out) < GWN(p_in). For an inverted face,
// GWN(p_out) > GWN(p_in). The DIFFERENCE is meaningful even when
// absolute GWNs are fractional (membrane interior).
//
// Algorithm details:
//   - Sample point: p = centroid + ε · face_normal, where ε is chosen
//     to be just larger than the local surface "thickness" (we use
//     bbox·1e-3 to clear local membrane / sliver pairs).
//   - The classification uses a robust threshold: flip if GWN > 0.5.
//     Robust against soup geometry where GWN can drift from {0, 1}.
//   - Patch-aware (orientation_patch BFS): faces in the same
//     consistent-winding patch flip together; the GWN is sampled at
//     ONE representative face per patch (the largest by area).
//     This reduces the GWN call count from O(F) to O(patches), making
//     the algorithm tractable on large meshes.
//
// Gating: only fires when the input has >= `min_antipar_pairs`
// anti-parallel manifold-edge pairs (the signal of the inverted-stripe
// pattern). For a mostly-correctly-oriented mesh the gate skips this
// stage at cost ~O(F) (one edge map build).
//
// Cost: O(P · F) where P is the patch count and F is the mesh size.
// For typical "inverted stripe" cases P is small (~2-5) so this is
// effectively O(F) per call.
//
// Caps:
//   - max_faces: if mesh.faces.size() > this, skip entirely. The O(P·F)
//     cost becomes prohibitive on >1 M F meshes; for those, the existing
//     destruction fallback path handles it.
//
// Guard: only applied if it strictly reduces the antipar-pair count.
// Otherwise the result is rolled back (the input mesh might already be
// optimally oriented; over-aggressive flipping would re-create
// inversions).
OrientWnResult orient_by_winding_number(const Mesh& mesh,
                                        uint32_t min_antipar_pairs = 50,
                                        size_t   max_faces = 1000000);

} // namespace meshseal::stages
