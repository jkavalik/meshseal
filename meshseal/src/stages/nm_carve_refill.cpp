#include "nm_carve_refill.h"
#include "holes.h"
#include "collapse_nm.h"
#include "../internal/diagnostics.h"
#include <climits>
#include <algorithm>
#include <array>
#include <cstdint>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace meshseal::stages {

namespace {

static uint64_t ekey(uint32_t a, uint32_t b) {
    if (a > b) std::swap(a, b);
    return (static_cast<uint64_t>(a) << 32) | b;
}

} // anonymous namespace

NmCarveRefillResult nm_carve_refill(const Mesh& mesh, int halo_rings) {
    NmCarveRefillResult result;
    result.mesh = mesh;
    if (mesh.faces.empty()) return result;

    // edge -> incident face list
    std::unordered_map<uint64_t, std::vector<uint32_t>> e2f;
    e2f.reserve(mesh.faces.size() * 3);
    for (uint32_t fi = 0; fi < mesh.faces.size(); ++fi) {
        const auto& f = mesh.faces[fi];
        for (int k = 0; k < 3; ++k) {
            e2f[ekey(f[k], f[(k+1)%3])].push_back(fi);
        }
    }
    // NM edges
    std::vector<uint64_t> nm_edges;
    for (const auto& kv : e2f)
        if (kv.second.size() > 2) nm_edges.push_back(kv.first);
    result.nm_before = static_cast<uint32_t>(nm_edges.size());
    if (nm_edges.empty()) return result;

    // Seed bad set with NM-incident faces.
    std::unordered_set<uint32_t> bad;
    for (uint64_t ek : nm_edges) {
        for (uint32_t fi : e2f[ek]) bad.insert(fi);
    }
    // vertex -> face map for halo expansion
    std::vector<std::vector<uint32_t>> v2f(mesh.vertices.size());
    for (uint32_t fi = 0; fi < mesh.faces.size(); ++fi) {
        const auto& f = mesh.faces[fi];
        for (int k = 0; k < 3; ++k) v2f[f[k]].push_back(fi);
    }
    for (int r = 0; r < halo_rings; ++r) {
        std::vector<uint32_t> add;
        for (uint32_t fi : bad) {
            const auto& f = mesh.faces[fi];
            for (int k = 0; k < 3; ++k)
                for (uint32_t nb : v2f[f[k]])
                    if (bad.find(nb) == bad.end()) add.push_back(nb);
        }
        for (uint32_t fi : add) bad.insert(fi);
    }
    result.faces_carved = static_cast<uint32_t>(bad.size());

    // Build carved mesh (drop bad faces; vertices kept as-is, unreferenced
    // are harmless and let downstream stages re-use stable indices).
    Mesh carved;
    carved.vertices = mesh.vertices;
    carved.faces.reserve(mesh.faces.size() - bad.size());
    for (uint32_t fi = 0; fi < mesh.faces.size(); ++fi) {
        if (bad.find(fi) == bad.end()) carved.faces.push_back(mesh.faces[fi]);
    }

    // Refill the holes opened by the carve via Liepa, then iterate
    // fill+collapse_nm until stable. Liepa often introduces new NM edges
    // when its fill triangles cross neighbouring geometry; collapse_nm
    // erases those. The combination converges in a few rounds.
    auto fr = stages::fill_holes(carved);
    result.faces_filled = fr.faces_added;
    Mesh patched = std::move(fr.mesh);
    int last_def = INT_MAX;
    for (int iter = 0; iter < 6; ++iter) {
        auto d = internal::compute_diagnostics(patched);
        const int cur_def = static_cast<int>(d.open_boundary_edges) +
                            static_cast<int>(d.non_manifold_edges);
        if (cur_def == 0 || cur_def >= last_def) break;
        last_def = cur_def;
        if (d.non_manifold_edges > 0) {
            auto cr = stages::collapse_nm_region(patched);
            if (cr.applied && cr.nm_after < cr.nm_before) {
                patched = std::move(cr.mesh);
            }
        }
        if (d.open_boundary_edges > 0) {
            auto fr2 = stages::fill_holes(patched);
            if (fr2.holes_filled > 0) {
                result.faces_filled += fr2.faces_added;
                patched = std::move(fr2.mesh);
            }
        }
    }

    result.mesh = std::move(patched);
    auto post = internal::compute_diagnostics(result.mesh);
    result.nm_after = static_cast<uint32_t>(post.non_manifold_edges);
    result.applied = true;
    return result;
}

} // namespace meshseal::stages
