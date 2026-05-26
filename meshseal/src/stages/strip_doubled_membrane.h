#pragma once
#include "../../include/meshseal/meshseal.h"
#include <cstdint>

namespace meshseal::stages {

struct StripDoubledMembraneResult {
    Mesh     mesh;
    bool     applied       = false;  // true iff the strip was adopted
    uint32_t pairs_found   = 0;      // antipar+proximity pairs detected
    uint32_t faces_removed = 0;      // = 2 * pairs_found when applied
    uint32_t total_faces   = 0;      // input face count
};

// Strip wholesale doubled membranes from the input.
//
// PROBLEM TARGETED. Multi-color / multi-shell exports from authoring tools
// (PrusaSlicer, Bambu Studio, MeshMixer) frequently contain an INTERNAL
// DOUBLED MEMBRANE: where two originally-separate solid bodies were left
// touching during export, the contact interface survives as two
// back-to-back face sheets with antiparallel normals. STL strips per-shell
// metadata, so the load-then-weld pipeline cannot distinguish "two solids
// that touched" from "one solid with a deliberate internal feature". The
// downstream `nm_edge`, `nm_vertex`, `orient` stages produce inverted-
// normal stripes + residual antipar manifold-edge pairs that look like
// real defects in the output.
//
// CRITERION. For each face F with centroid C_F and unit normal N_F, find
// any face F' such that
//     ||C_F - C_F'|| < tol_mul * bbox_diag
//     N_F · N_F'     < antipar_dot_max
// (default tol_mul = 2e-3, antipar_dot_max = -0.85.) Mark both as
// "membrane" if neither is already marked (greedy pairing).
//
// GATE. The detected fraction `2 * pairs_found / total_faces` must be at
// least `min_fraction` (default 0.10 = 10%) for the strip to be adopted.
// This separates "wholesale doubled-membrane" inputs (Bee_v3 STL: 81%
// detected; vol[3] alone: ~24% — both trigger) from "small localised
// doubled patch" inputs (trumpet: 3% detected — does NOT trigger, leaves
// the existing collapse_nm pipeline to handle it). The gate is what makes
// the stage safe to run unconditionally: clean inputs detect 0% and
// pass through untouched.
//
// VALIDATED ON Bee_v3 (vcela/Bee_v3.stl):
//   - Whole STL:  F=137588, nm=82924, antipar=31  → meshseal baseline
//                  F=78920, nm=44, antipar=1074 (and c=2, fragmented)
//                 → with strip: F≈23000, nm=0, antipar=22, c=1, vol matched
//                  to sum of per-shell volumes.
//   - Per-volume vol[3] alone: input antipar=10 → baseline output antipar=261
//                              → with strip: antipar=5, vol preserved.
//
// The stage is structural, not heuristic: a face is dropped only when it
// has a back-to-back antiparallel partner within ε. There is no smoothing,
// no remeshing, no winding-number reasoning — just position+normal-based
// pair identification + removal.
// `min_total_faces` guards tiny meshes (a 2-face back-to-back input or
// the 24-face big-and-tiny-cube test fixture would otherwise read as 100%
// or 50% membrane and falsely trigger the strip on legitimate small
// objects). The default of 100 lets all corpus fixtures except synthetic
// degenerate probes pass the gate.
//
// `max_fraction` guards FULLY-doubled-soup inputs: scan reconstructions
// (e.g. room_surface_mesh) have a back-to-back partner for *every* face
// (because the scanner picks up both sides of every surface) -> 100% of
// faces detected as pairs -> strip would erase the entire mesh.
// Such inputs are not "internal membrane" cases; they need the existing
// soup_reconstruct pathway. The default 0.90 lets Bee_v3 (81 %) through
// while skipping fully-doubled soups.
StripDoubledMembraneResult strip_doubled_membrane(
    const Mesh& mesh,
    double   tol_mul         = 2e-3,
    double   antipar_dot_max = -0.85,
    double   min_fraction    = 0.10,
    double   max_fraction    = 0.90,
    uint32_t min_total_faces = 100);

} // namespace meshseal::stages
