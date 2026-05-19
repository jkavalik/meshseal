#pragma once
#include "../../include/meshseal/meshseal.h"
#include <cstdint>

namespace meshseal::stages {

struct TJunctionResult {
    Mesh     mesh;
    uint32_t edges_split = 0;   // T-junction edges resolved
    uint32_t faces_split = 0;   // faces fanned through collinear vertices
};

// Resolve T-junctions.
//
// A T-junction is a non-conforming edge: an edge (a,b) carried by one face
// while the other side of the surface is split — one or more mesh vertices
// lie collinear on the *interior* of segment (a,b), and the faces there use
// the sub-edges (a,v1), (v1,v2), …, (vk,b) instead. The spanning edge has no
// 2-to-1 partner, so it reads as a boundary (or non-manifold) edge and
// hole-fill cannot close it (a zero-width collinear slit is a zero-area
// loop; Liepa rejects it).
//
// Fix: for every problem edge (count != 2) that has collinear interior
// vertices, fan-split each incident face through those vertices — face
// (a,b,c) → (a,v1,c),(v1,v2,c),…,(vk,b,c). The sub-edges then pair 2-to-1
// with the other side and the slit closes with no degenerate faces.
//
// Resolves real-world collinear slits (e.g. ore's bnd3 residual) and the
// `fill_bad_diagonal` fixture (a cube with a non-conforming top diagonal).
TJunctionResult split_tjunctions(const Mesh& mesh);

} // namespace meshseal::stages
