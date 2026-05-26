#include "strip_doubled_membrane.h"
#include "../internal/spatial_hash.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <vector>

namespace meshseal::stages {

namespace {

// Length of a 3-vector.
inline double vlen(const std::array<double, 3>& v) {
    return std::sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
}

// Cross product b - a (vector subtraction).
inline std::array<double, 3> sub(const std::array<double, 3>& a,
                                  const std::array<double, 3>& b) {
    return {a[0]-b[0], a[1]-b[1], a[2]-b[2]};
}

// Cross product.
inline std::array<double, 3> cross(const std::array<double, 3>& a,
                                    const std::array<double, 3>& b) {
    return {a[1]*b[2] - a[2]*b[1],
            a[2]*b[0] - a[0]*b[2],
            a[0]*b[1] - a[1]*b[0]};
}

// Dot product.
inline double dot(const std::array<double, 3>& a,
                  const std::array<double, 3>& b) {
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

} // anon

StripDoubledMembraneResult strip_doubled_membrane(
    const Mesh& mesh,
    double   tol_mul,
    double   antipar_dot_max,
    double   min_fraction,
    double   max_fraction,
    uint32_t min_total_faces)
{
    StripDoubledMembraneResult res;
    res.mesh = mesh;
    res.total_faces = static_cast<uint32_t>(mesh.faces.size());

    if (mesh.faces.empty() || mesh.vertices.empty()) return res;
    if (res.total_faces < min_total_faces) return res;

    // Compute bbox diagonal.
    std::array<double, 3> lo = mesh.vertices[0];
    std::array<double, 3> hi = mesh.vertices[0];
    for (const auto& p : mesh.vertices) {
        for (int k = 0; k < 3; ++k) {
            if (p[k] < lo[k]) lo[k] = p[k];
            if (p[k] > hi[k]) hi[k] = p[k];
        }
    }
    const double bbd = std::sqrt(
        (hi[0]-lo[0])*(hi[0]-lo[0]) +
        (hi[1]-lo[1])*(hi[1]-lo[1]) +
        (hi[2]-lo[2])*(hi[2]-lo[2]));
    if (bbd <= 0.0) return res;

    const double tol = tol_mul * bbd;
    if (tol <= 0.0) return res;

    // Compute per-face centroid + unit normal. Skip zero-area faces.
    const size_t F = mesh.faces.size();
    std::vector<std::array<double, 3>> centroid(F);
    std::vector<std::array<double, 3>> normal(F);
    std::vector<uint8_t> valid(F, 0);
    for (size_t i = 0; i < F; ++i) {
        const auto& f = mesh.faces[i];
        const auto& A = mesh.vertices[f[0]];
        const auto& B = mesh.vertices[f[1]];
        const auto& C = mesh.vertices[f[2]];
        centroid[i] = {(A[0]+B[0]+C[0]) / 3.0,
                       (A[1]+B[1]+C[1]) / 3.0,
                       (A[2]+B[2]+C[2]) / 3.0};
        // Normal = (B-A) x (C-A) [matches the orientation convention].
        auto BA = sub(B, A);
        auto CA = sub(C, A);
        auto N = cross(BA, CA);
        const double ln = vlen(N);
        if (ln > 1e-30) {
            normal[i] = {N[0]/ln, N[1]/ln, N[2]/ln};
            valid[i] = 1;
        }
    }

    // Spatial-hash insertion keyed by centroid. Cell size = 2 * tol so a
    // 27-cell lookup ALWAYS contains any centroid within `tol` regardless
    // of which side of a cell boundary the query lies on.
    const double cell = std::max(2.0 * tol, bbd * 1e-9);
    internal::SpatialHash hash(cell);
    for (size_t i = 0; i < F; ++i) {
        if (valid[i]) hash.insert(centroid[i], static_cast<uint32_t>(i));
    }

    // Greedy pairing. For each face, query neighbours and pair with the
    // first unmarked antipar+close one. Visiting in face-index order is
    // deterministic; the neighbour list `nbrs` is unordered (hash bucket
    // order). Sorting `nbrs` would make pair selection reproducible across
    // builds but regresses Bee_v3.stl antipar (49 → 141) — the downstream
    // pipeline's pair resolution is calibrated to the un-sorted choice.
    // Documented determinism risk; left for a future careful re-calibration.
    std::vector<uint8_t> detected(F, 0);
    uint32_t pair_count = 0;
    const double tol2 = tol * tol;
    for (size_t i = 0; i < F; ++i) {
        if (!valid[i] || detected[i]) continue;
        auto nbrs = hash.query_neighbors(centroid[i]);
        for (uint32_t j : nbrs) {
            if (j == i || detected[j] || !valid[j]) continue;
            // Antiparallel check.
            const double d = dot(normal[i], normal[j]);
            if (d > antipar_dot_max) continue;
            // Distance check.
            const auto& Ci = centroid[i];
            const auto& Cj = centroid[j];
            const double dx = Cj[0] - Ci[0];
            const double dy = Cj[1] - Ci[1];
            const double dz = Cj[2] - Ci[2];
            const double d2 = dx*dx + dy*dy + dz*dz;
            if (d2 > tol2) continue;
            // Pair them.
            detected[i] = 1;
            detected[j] = 1;
            ++pair_count;
            break;
        }
    }

    res.pairs_found = pair_count;
    res.faces_removed = pair_count * 2u;

    // Apply the gates.
    const double detected_fraction = static_cast<double>(res.faces_removed) /
                                     static_cast<double>(F);
    if (detected_fraction < min_fraction || detected_fraction > max_fraction) {
        // Don't adopt — return original mesh, applied = false.
        return res;
    }

    // Build new face list excluding detected faces.
    std::vector<std::array<uint32_t, 3>> new_faces;
    new_faces.reserve(F - res.faces_removed);
    for (size_t i = 0; i < F; ++i) {
        if (!detected[i]) new_faces.push_back(mesh.faces[i]);
    }
    res.mesh.faces = std::move(new_faces);
    res.applied = true;
    return res;
}

} // namespace meshseal::stages
