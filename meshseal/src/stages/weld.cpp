#include "weld.h"
#include "../internal/spatial_hash.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

namespace meshseal::stages {

WeldResult weld(const Mesh& mesh, std::optional<double> tolerance) {
    const uint32_t n_verts = static_cast<uint32_t>(mesh.vertices.size());
    const uint32_t n_faces = static_cast<uint32_t>(mesh.faces.size());

    // ------------------------------------------------------------------ //
    // 1. Compute tolerance (if not provided)
    // ------------------------------------------------------------------ //
    double tol;
    if (tolerance.has_value()) {
        tol = tolerance.value();
    } else {
        // Bounding box diagonal
        double lo[3] = { std::numeric_limits<double>::max(),
                         std::numeric_limits<double>::max(),
                         std::numeric_limits<double>::max() };
        double hi[3] = { std::numeric_limits<double>::lowest(),
                         std::numeric_limits<double>::lowest(),
                         std::numeric_limits<double>::lowest() };
        for (const auto& v : mesh.vertices) {
            for (int k = 0; k < 3; ++k) {
                lo[k] = std::min(lo[k], v[k]);
                hi[k] = std::max(hi[k], v[k]);
            }
        }
        double d = 0.0;
        if (n_verts > 0) {
            double dx = hi[0] - lo[0];
            double dy = hi[1] - lo[1];
            double dz = hi[2] - lo[2];
            d = std::sqrt(dx * dx + dy * dy + dz * dz);
        }

        // Mean edge length from first min(1000, F) faces
        double edge_sum = 0.0;
        uint32_t edge_count = 0;
        const uint32_t sample_faces = std::min(n_faces, static_cast<uint32_t>(1000));
        for (uint32_t f = 0; f < sample_faces; ++f) {
            const auto& face = mesh.faces[f];
            for (int e = 0; e < 3; ++e) {
                const auto& a = mesh.vertices[face[e]];
                const auto& b = mesh.vertices[face[(e + 1) % 3]];
                double ex = b[0] - a[0];
                double ey = b[1] - a[1];
                double ez = b[2] - a[2];
                edge_sum += std::sqrt(ex * ex + ey * ey + ez * ez);
                ++edge_count;
            }
        }
        const double mean_edge = (edge_count > 0) ? (edge_sum / edge_count) : 0.0;

        tol = std::max({ 1e-9, mean_edge * 1e-4, d * 1e-7 });
    }

    // ------------------------------------------------------------------ //
    // 2. Build SpatialHash and insert all vertices
    // ------------------------------------------------------------------ //
    internal::SpatialHash hash(tol);
    for (uint32_t i = 0; i < n_verts; ++i) {
        hash.insert(mesh.vertices[i], i);
    }

    // ------------------------------------------------------------------ //
    // 3. Build canonical vertex map (old_idx → canonical_idx)
    // ------------------------------------------------------------------ //
    // canonical_of[i] = the representative index for vertex i
    std::vector<uint32_t> canonical_of(n_verts);
    for (uint32_t i = 0; i < n_verts; ++i) canonical_of[i] = i;
    // new_index[i] = index in output vertex list for the representative i
    // Use n_verts as sentinel for "not yet assigned"
    std::vector<uint32_t> new_index(n_verts, n_verts);

    std::vector<std::array<double, 3>> out_vertices;
    out_vertices.reserve(n_verts);

    for (uint32_t i = 0; i < n_verts; ++i) {
        const auto& pi = mesh.vertices[i];
        const auto neighbors = hash.query_neighbors(pi);

        // Find the smallest-indexed neighbor within tolerance
        uint32_t best = i;
        for (uint32_t nb : neighbors) {
            if (nb < best) {
                const auto& pnb = mesh.vertices[nb];
                double dx = pnb[0] - pi[0];
                double dy = pnb[1] - pi[1];
                double dz = pnb[2] - pi[2];
                if (std::sqrt(dx * dx + dy * dy + dz * dz) <= tol) {
                    best = nb;
                }
            }
        }
        // Resolve chain: walk to the true root (safe — best < i, already processed)
        while (canonical_of[best] != best) best = canonical_of[best];
        canonical_of[i] = best;

        // If i itself is the representative, assign a new output index
        if (best == i) {
            new_index[i] = static_cast<uint32_t>(out_vertices.size());
            out_vertices.push_back(pi);
        }
    }

    // ------------------------------------------------------------------ //
    // 4. Remap faces, dropping degenerate ones
    // ------------------------------------------------------------------ //
    std::vector<std::array<uint32_t, 3>> out_faces;
    out_faces.reserve(n_faces);

    for (const auto& face : mesh.faces) {
        const uint32_t a = new_index[canonical_of[face[0]]];
        const uint32_t b = new_index[canonical_of[face[1]]];
        const uint32_t c = new_index[canonical_of[face[2]]];

        if (a == b || b == c || a == c) {
            continue; // degenerate after welding
        }
        out_faces.push_back({ a, b, c });
    }

    // ------------------------------------------------------------------ //
    // 5. Return result
    // ------------------------------------------------------------------ //
    WeldResult result;
    result.mesh.vertices = std::move(out_vertices);
    result.mesh.faces    = std::move(out_faces);
    result.vertices_before = n_verts;
    result.vertices_after  = static_cast<uint32_t>(result.mesh.vertices.size());
    result.tolerance_used  = tol;
    return result;
}

} // namespace meshseal::stages

