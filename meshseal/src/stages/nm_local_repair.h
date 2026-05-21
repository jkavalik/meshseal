#pragma once
#include "../../include/meshseal/meshseal.h"
#include <cstdint>

namespace meshseal::stages {

struct NmLocalRepairResult {
    Mesh     mesh;
    uint32_t merges        = 0;  // NM-local proximity welds applied
    uint32_t pairs_removed = 0;  // exact back-to-back face pairs removed
    uint32_t faces_dropped = 0;  // total face count drop (degenerates + pairs)
};

// Two combined operations targeted at residual non-manifold edges left
// over after collapse_nm:
//
//   1. **NM-local proximity weld.** For each NM-edge endpoint, merge any
//      1-ring neighbour that sits within `tol_rel * bbox_diag` of it.
//      Catches design-defect near-coincident duplicate vertices that the
//      global metric weld (typically bbox*1e-7) is too tight to handle —
//      e.g. t10k_1582375's chamfer corners with vertices 0.004 apart that
//      should be one vertex.
//
//   2. **Strict back-to-back face dedup.** After (1), find face pairs
//      with IDENTICAL vertex sets (modulo cyclic permutation) and
//      OPPOSITE winding. Remove both. This is what remove_thin_features
//      tries to do with a centroid+normal metric, but the metric also
//      matches nearby non-pair faces and worsens the result on inputs
//      like bumpy_white. The strict vertex-set match catches only true
//      duplicates.
//
// Combined into one stage because (1) often creates exact back-to-back
// pairs that (2) must clean up afterwards. The whole stage is rolled
// back if the result does not strictly improve (bnd+nm) or moves
// signed volume by more than 0.5 % — same envelope as collapse_nm.
NmLocalRepairResult nm_local_repair(const Mesh& mesh,
                                    double tol_rel = 1e-4);

} // namespace meshseal::stages
