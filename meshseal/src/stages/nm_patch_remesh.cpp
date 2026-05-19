#include "nm_patch_remesh.h"
#include "holes.h"
#include "../internal/diagnostics.h"
#include <array>
#include <cmath>
#include <cstdint>
#include <unordered_map>
#include <vector>

namespace meshseal::stages {

namespace {

uint64_t ekey(uint32_t a, uint32_t b) {
    if (a > b) std::swap(a, b);
    return (static_cast<uint64_t>(a) << 32) | b;
}

// Unit normal of face fi (zero vector if degenerate).
std::array<double, 3> face_normal(const Mesh& mesh, uint32_t fi) {
    const auto& f = mesh.faces[fi];
    const auto& a = mesh.vertices[f[0]];
    const auto& b = mesh.vertices[f[1]];
    const auto& c = mesh.vertices[f[2]];
    double ux = b[0]-a[0], uy = b[1]-a[1], uz = b[2]-a[2];
    double wx = c[0]-a[0], wy = c[1]-a[1], wz = c[2]-a[2];
    double nx = uy*wz - uz*wy, ny = uz*wx - ux*wz, nz = ux*wy - uy*wx;
    double L = std::sqrt(nx*nx + ny*ny + nz*nz);
    if (L < 1e-300) return {0.0, 0.0, 0.0};
    return {nx/L, ny/L, nz/L};
}

} // anonymous namespace

NmPatchResult remesh_nm_patches(const Mesh& mesh, int rings,
                                bool coplanar_expand, double max_patch_frac) {
    NmPatchResult result;
    result.mesh = mesh;

    const uint32_t nf = static_cast<uint32_t>(mesh.faces.size());
    if (nf == 0) return result;

    // Edge → incident faces.
    std::unordered_map<uint64_t, std::vector<uint32_t>> e2f;
    e2f.reserve(nf * 3);
    for (uint32_t fi = 0; fi < nf; ++fi) {
        const auto& f = mesh.faces[fi];
        e2f[ekey(f[0], f[1])].push_back(fi);
        e2f[ekey(f[1], f[2])].push_back(fi);
        e2f[ekey(f[2], f[0])].push_back(fi);
    }

    // Seed the patch with every face incident to a non-manifold edge.
    std::vector<bool> in_patch(nf, false);
    uint32_t nm_edges = 0;
    for (const auto& kv : e2f) {
        if (kv.second.size() > 2) {
            ++nm_edges;
            for (uint32_t fi : kv.second) in_patch[fi] = true;
        }
    }
    result.nm_before = nm_edges;
    if (nm_edges == 0) return result;  // nothing to do; applied stays false
    result.applied = true;

    // Coplanar expansion: flood-fill each seed face into its whole coplanar
    // edge-connected sheet. |n_i·n_j| > 0.999 — anti-parallel counts too, so
    // both layers of a doubled flat membrane are swallowed. Capped so a
    // runaway sheet (entire flat wall) cannot consume the model.
    if (coplanar_expand) {
        const uint32_t cap = static_cast<uint32_t>(max_patch_frac * nf);
        std::vector<std::array<double,3>> nrm(nf);
        for (uint32_t fi = 0; fi < nf; ++fi) nrm[fi] = face_normal(mesh, fi);

        std::vector<uint32_t> stack;
        for (uint32_t fi = 0; fi < nf; ++fi)
            if (in_patch[fi]) stack.push_back(fi);
        uint32_t patch_count = static_cast<uint32_t>(stack.size());

        while (!stack.empty() && patch_count < cap) {
            uint32_t fi = stack.back();
            stack.pop_back();
            const auto& ni = nrm[fi];
            if (ni[0]==0.0 && ni[1]==0.0 && ni[2]==0.0) continue;
            const auto& f = mesh.faces[fi];
            for (int k = 0; k < 3; ++k) {
                for (uint32_t nb : e2f[ekey(f[k], f[(k+1)%3])]) {
                    if (in_patch[nb]) continue;
                    const auto& nn = nrm[nb];
                    double d = ni[0]*nn[0] + ni[1]*nn[1] + ni[2]*nn[2];
                    if (std::abs(d) > 0.999) {
                        in_patch[nb] = true;
                        ++patch_count;
                        stack.push_back(nb);
                        if (patch_count >= cap) break;
                    }
                }
                if (patch_count >= cap) break;
            }
        }
    }

    // BFS-expand the patch by `rings` edge-connected steps.
    for (int r = 0; r < rings; ++r) {
        std::vector<uint32_t> frontier;
        for (uint32_t fi = 0; fi < nf; ++fi)
            if (in_patch[fi]) frontier.push_back(fi);
        for (uint32_t fi : frontier) {
            const auto& f = mesh.faces[fi];
            for (int k = 0; k < 3; ++k) {
                for (uint32_t nb : e2f[ekey(f[k], f[(k+1)%3])])
                    in_patch[nb] = true;
            }
        }
    }

    // Build the patch-deleted mesh (keep all vertices — unused ones are
    // harmless and keeping indices stable lets fill_holes append cleanly).
    Mesh holed;
    holed.vertices = mesh.vertices;
    holed.faces.reserve(nf);
    for (uint32_t fi = 0; fi < nf; ++fi) {
        if (in_patch[fi]) { ++result.patch_removed; continue; }
        holed.faces.push_back(mesh.faces[fi]);
    }

    // Close the opened boundary loop(s) with Liepa hole-fill.
    const uint32_t faces_before_fill = static_cast<uint32_t>(holed.faces.size());
    auto hr = stages::fill_holes(holed);
    result.faces_filled = static_cast<uint32_t>(hr.mesh.faces.size()) - faces_before_fill;

    auto diag = internal::compute_diagnostics(hr.mesh);
    result.nm_after = static_cast<uint32_t>(diag.non_manifold_edges);
    result.mesh = std::move(hr.mesh);
    return result;
}

} // namespace meshseal::stages
