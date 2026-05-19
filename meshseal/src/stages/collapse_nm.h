#pragma once
#include "../../include/meshseal/meshseal.h"
#include <cstdint>

namespace meshseal::stages {

struct CollapseNmResult {
    Mesh     mesh;
    bool     applied    = false;  // false if there were no NM edges
    uint32_t nm_before  = 0;
    uint32_t nm_after   = 0;
    uint32_t collapses  = 0;      // number of edge collapses performed
};

// Progressive guarded edge-collapse to erase non-manifold flap regions.
//
// The residual non-manifold geometry on real-world fixtures (trumpet, ...) is
// a degenerate near-zero-volume doubled membrane — a CSG-union seam where two
// near-coplanar sheets were left stacked. Deleting + hole-filling it fails
// (the refill re-creates the NM). Instead, this stage *erases* the membrane
// by progressively collapsing its edges:
//
//   repeat: pick the shortest edge in the NM region and collapse it, where
//   every collapse is GUARDED so the mesh never gets worse —
//     * non-manifold-edge count must not rise,
//     * open-boundary-edge count must not rise (watertight preserved at every
//       step — the only way boundary grows is a genuine rip, which is barred),
//     * no real (non-sliver) face may have its normal flipped.
//   Two tiers: an "improve" collapse strictly drops the NM count; a "neutral"
//   collapse keeps NM/boundary equal while shrinking the region (decimating
//   the degenerate patch to escape a stall). Stop when NM==0 or no guarded
//   collapse exists.
//
// Because boundary can never grow, the operation only ever erases an internal
// degenerate pocket; geometry change is tiny and local (empirically <0.01 %
// volume on trumpet). `max_collapse_len` (auto = 0.05·bbox_diag when <0)
// caps which edges may be collapsed so real feature edges are never touched.
CollapseNmResult collapse_nm_region(const Mesh& mesh,
                                    double max_collapse_len = -1.0,
                                    int    max_collapses    = 4000);

} // namespace meshseal::stages
