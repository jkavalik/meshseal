#pragma once
#include "../../include/meshseal/meshseal.h"
#include <cstdint>

namespace meshseal::stages {

struct NmCarveRefillResult {
    Mesh     mesh;
    uint32_t faces_carved = 0;
    uint32_t faces_filled = 0;
    uint32_t nm_before    = 0;
    uint32_t nm_after     = 0;
    bool     applied      = false;
};

// Resolve stubborn residual non-manifold edges by carving the local
// region around each NM edge and refilling with Liepa hole-fill.
//
// For each NM edge in the mesh: collect all incident faces, then expand
// by `halo_rings` rings of vertex-neighbour faces (default 1 — minimal
// halo gives the smallest acceptable carve in practice). Delete those
// faces. The result is a mesh with small open boundaries where the NM
// region used to be. Caller (or this stage) then runs `fill_holes` to
// close them, producing clean topology.
//
// Tested on real-world scan inputs where NM edges sit inside a cluster
// of self-intersections (black_vase): the local carve + refill resolves
// the NM without touching the rest of the mesh. black_vase: 41 faces
// carved, vol +0.03%, nm 1 -> 0.
//
// NOTE: this is *not* an SI cleanup stage on its own. It targets the
// case where collapse_nm, nm_patch_remesh and nm_local_repair have all
// failed to resolve a stubborn residual nm. Apply LATE in the pipeline.
NmCarveRefillResult nm_carve_refill(const Mesh& mesh, int halo_rings = 1);

} // namespace meshseal::stages
