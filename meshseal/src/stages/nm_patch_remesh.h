#pragma once
#include "../../include/meshseal/meshseal.h"
#include <cstdint>

namespace meshseal::stages {

struct NmPatchResult {
    Mesh     mesh;
    bool     applied        = false;  // false if there were no NM edges
    uint32_t nm_before      = 0;
    uint32_t nm_after       = 0;
    uint32_t patch_removed  = 0;      // faces deleted from the NM patch
    uint32_t faces_filled   = 0;      // faces added by the refill
};

// Local remesh of non-manifold patches (Idea C).
//
// A non-manifold edge with ≥3 incident faces cannot be made 2-manifold by
// edge tricks when the cause is structural — a doubled-surface flap (bad CSG
// self-intersection resolution) or a multi-sheet pinch. The robust fix is
// LOCAL: delete the small patch of faces around each NM edge, then close the
// resulting boundary loop(s) with Liepa hole-fill, which produces a clean
// single-layer 2-manifold disk.
//
// Algorithm:
//   1. find NM edges (≥3 incident faces);
//   2. seed the patch with all NM-incident faces;
//   3. if `coplanar_expand`, flood-fill each seed face into its full
//      coplanar edge-connected sheet (|n_i·n_j| > ~0.999) — this swallows a
//      whole doubled flat membrane (e.g. a trumpet bell-rim disc that is
//      modelled as two coincident independently-triangulated layers, which a
//      tiny BFS patch cannot resolve because Liepa just re-creates the NM);
//   4. BFS-expand `rings` further edge-connected steps;
//   5. delete the patch faces — opens boundary loop(s);
//   6. `fill_holes` (Liepa) re-triangulates the loop(s) as clean disks.
//
// `rings` controls extra patch growth beyond the coplanar sheet. The caller
// should apply a total-defect keep/reject guard. `max_patch_frac` caps the
// coplanar flood-fill so a runaway sheet (whole flat wall) is not swallowed.
NmPatchResult remesh_nm_patches(const Mesh& mesh,
                                int    rings           = 1,
                                bool   coplanar_expand = false,
                                double max_patch_frac  = 0.40);

} // namespace meshseal::stages
