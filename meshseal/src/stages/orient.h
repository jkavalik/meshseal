#pragma once
#include "../../include/meshseal/meshseal.h"
#include <cstdint>
#include <vector>

namespace meshseal::stages {

struct OrientResult {
    Mesh mesh;
    uint32_t faces_flipped;
    uint32_t components_flipped; // components where all faces had winding reversed
    bool was_orientable;         // false if mesh has non-orientable regions
};

// Repair face winding order so all normals point outward.
// Algorithm:
//   1. BFS per connected component: propagate consistent winding
//   2. For each component, compute signed volume
//   3. If signed volume < 0, flip all faces in that component
//
// `do_signed_volume_flip` controls step 3. Pass `false` when re-orienting
// after a topology-changing stage (e.g. nm_edge fin removal): the BFS pass
// is safe to re-run — it just propagates winding through newly-manifold
// seams — but the signed-volume flip can wrongly invert small components
// whose sign happened to flip because of face removal.
OrientResult orient_mesh(const Mesh& mesh, bool do_signed_volume_flip = true);

} // namespace meshseal::stages
