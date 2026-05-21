#include "si_split.h"
#include "../internal/arrangement.h"
#include "../internal/vec3.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace meshseal::stages {

namespace {

using internal::Vec3d;

static uint64_t cell_key(int32_t cx, int32_t cy, int32_t cz) {
    auto pack = [](int32_t v) -> uint64_t {
        return static_cast<uint64_t>(static_cast<int64_t>(v) + (1 << 20)) & 0x1FFFFF;
    };
    return pack(cx) | (pack(cy) << 21) | (pack(cz) << 42);
}

// Cut consolidation: Möller returns the intersection segment clipped to
// BOTH triangles of a pair, so the same physical chord through a triangle
// may be reported by several pairs as partial sub-segments. Merge cuts on
// the same triangle that share a line, within tolerance.
static std::vector<std::pair<Vec3d, Vec3d>> consolidate(
    const std::vector<std::pair<Vec3d, Vec3d>>& cuts,
    double line_tol)
{
    if (cuts.size() <= 1) return cuts;
    const size_t n = cuts.size();
    std::vector<std::pair<Vec3d, Vec3d>> merged;
    std::vector<bool> taken(n, false);
    for (size_t i = 0; i < n; ++i) {
        if (taken[i]) continue;
        const Vec3d& A   = cuts[i].first;
        const Vec3d  AB  = internal::sub(cuts[i].second, A);
        const double L   = internal::norm(AB);
        if (L < 1e-12) { taken[i] = true; continue; }
        const Vec3d  dir = internal::scale(AB, 1.0 / L);
        double t_lo = 0.0, t_hi = L;
        Vec3d  pt_lo = A;
        Vec3d  pt_hi = cuts[i].second;
        taken[i] = true;
        for (size_t j = i + 1; j < n; ++j) {
            if (taken[j]) continue;
            auto on_line = [&](const Vec3d& P) {
                const Vec3d AP   = internal::sub(P, A);
                const double t   = internal::dot(AP, dir);
                const Vec3d proj = internal::add(A, internal::scale(dir, t));
                const Vec3d d    = internal::sub(P, proj);
                return internal::norm(d) < line_tol;
            };
            if (!on_line(cuts[j].first) || !on_line(cuts[j].second)) continue;
            const double tA = internal::dot(internal::sub(cuts[j].first, A), dir);
            const double tB = internal::dot(internal::sub(cuts[j].second, A), dir);
            const double tmin = std::min(tA, tB);
            const double tmax = std::max(tA, tB);
            if (tmin < t_lo) { t_lo = tmin; pt_lo = (tA < tB) ? cuts[j].first : cuts[j].second; }
            if (tmax > t_hi) { t_hi = tmax; pt_hi = (tA > tB) ? cuts[j].first : cuts[j].second; }
            taken[j] = true;
        }
        merged.emplace_back(pt_lo, pt_hi);
    }
    return merged;
}

} // anonymous namespace

SiSplitResult resolve_self_intersections(const Mesh& mesh) {
    SiSplitResult result;
    result.faces_before = static_cast<uint32_t>(mesh.faces.size());
    result.faces_after  = result.faces_before;
    result.mesh         = mesh;
    if (mesh.faces.size() < 2 || mesh.vertices.empty()) return result;

    const auto& V = mesh.vertices;
    const auto& F = mesh.faces;
    const size_t nf = F.size();

    // bbox + median edge → grid cell
    Vec3d lo = V[0], hi = V[0];
    for (const auto& p : V)
        for (int k = 0; k < 3; ++k) {
            lo[k] = std::min(lo[k], p[k]);
            hi[k] = std::max(hi[k], p[k]);
        }
    const double dx = hi[0]-lo[0], dy = hi[1]-lo[1], dz = hi[2]-lo[2];
    const double bbox_diag = std::sqrt(dx*dx + dy*dy + dz*dz);
    if (bbox_diag < 1e-30) return result;
    std::vector<double> edge_lens;
    edge_lens.reserve(nf * 3);
    for (const auto& t : F) {
        for (int k = 0; k < 3; ++k) {
            const Vec3d d = internal::sub(V[t[k]], V[t[(k+1)%3]]);
            edge_lens.push_back(internal::norm(d));
        }
    }
    std::nth_element(edge_lens.begin(),
                     edge_lens.begin() + edge_lens.size()/2,
                     edge_lens.end());
    double cell = edge_lens[edge_lens.size()/2];
    if (cell < bbox_diag * 1e-6) cell = bbox_diag * 1e-6;

    // Per-face triangle bbox + bucket grid.
    std::vector<std::array<Vec3d, 2>> tri_bb(nf);
    std::unordered_map<uint64_t, std::vector<uint32_t>> grid;
    grid.reserve(nf * 2);
    for (size_t fi = 0; fi < nf; ++fi) {
        const auto& t = F[fi];
        Vec3d mn = V[t[0]], mx = mn;
        for (int k = 1; k < 3; ++k)
            for (int j = 0; j < 3; ++j) {
                mn[j] = std::min(mn[j], V[t[k]][j]);
                mx[j] = std::max(mx[j], V[t[k]][j]);
            }
        tri_bb[fi] = { mn, mx };
        const int32_t cx0 = static_cast<int32_t>(std::floor(mn[0] / cell));
        const int32_t cy0 = static_cast<int32_t>(std::floor(mn[1] / cell));
        const int32_t cz0 = static_cast<int32_t>(std::floor(mn[2] / cell));
        const int32_t cx1 = static_cast<int32_t>(std::floor(mx[0] / cell));
        const int32_t cy1 = static_cast<int32_t>(std::floor(mx[1] / cell));
        const int32_t cz1 = static_cast<int32_t>(std::floor(mx[2] / cell));
        for (int32_t cx = cx0; cx <= cx1; ++cx)
        for (int32_t cy = cy0; cy <= cy1; ++cy)
        for (int32_t cz = cz0; cz <= cz1; ++cz)
            grid[cell_key(cx, cy, cz)].push_back(static_cast<uint32_t>(fi));
    }

    auto shares_vertex = [&](uint32_t i, uint32_t j) {
        const auto& ti = F[i]; const auto& tj = F[j];
        for (int a = 0; a < 3; ++a)
            for (int b = 0; b < 3; ++b)
                if (ti[a] == tj[b]) return true;
        return false;
    };
    auto bb_overlap = [](const std::array<Vec3d,2>& a, const std::array<Vec3d,2>& b) {
        for (int k = 0; k < 3; ++k)
            if (a[1][k] < b[0][k] || b[1][k] < a[0][k]) return false;
        return true;
    };

    // Collect per-face cut segments.
    std::vector<std::vector<std::pair<Vec3d, Vec3d>>> tri_cuts(nf);
    std::unordered_set<uint64_t> seen_pairs;
    seen_pairs.reserve(nf);
    for (size_t fi = 0; fi < nf; ++fi) {
        const auto& bb = tri_bb[fi];
        const int32_t cx0 = static_cast<int32_t>(std::floor(bb[0][0] / cell));
        const int32_t cy0 = static_cast<int32_t>(std::floor(bb[0][1] / cell));
        const int32_t cz0 = static_cast<int32_t>(std::floor(bb[0][2] / cell));
        const int32_t cx1 = static_cast<int32_t>(std::floor(bb[1][0] / cell));
        const int32_t cy1 = static_cast<int32_t>(std::floor(bb[1][1] / cell));
        const int32_t cz1 = static_cast<int32_t>(std::floor(bb[1][2] / cell));
        for (int32_t cx = cx0; cx <= cx1; ++cx)
        for (int32_t cy = cy0; cy <= cy1; ++cy)
        for (int32_t cz = cz0; cz <= cz1; ++cz) {
            auto it = grid.find(cell_key(cx, cy, cz));
            if (it == grid.end()) continue;
            for (uint32_t fj : it->second) {
                if (fj <= fi) continue;
                uint64_t pk = (static_cast<uint64_t>(fi) << 32) | fj;
                if (!seen_pairs.insert(pk).second) continue;
                if (!bb_overlap(tri_bb[fi], tri_bb[fj])) continue;
                if (shares_vertex(static_cast<uint32_t>(fi), fj)) continue;
                const auto& ti = F[fi]; const auto& tj = F[fj];
                auto seg = internal::tri_tri_intersect(
                    V[ti[0]], V[ti[1]], V[ti[2]],
                    V[tj[0]], V[tj[1]], V[tj[2]]);
                if (seg.size() == 2) {
                    tri_cuts[fi].emplace_back(seg[0], seg[1]);
                    tri_cuts[fj].emplace_back(seg[0], seg[1]);
                    ++result.pairs_found;
                }
            }
        }
    }
    if (result.pairs_found == 0) return result;

    // Consolidate co-linear cuts per face, then split each face.
    const double line_tol = cell * 0.02;
    const double edge_tol = cell * 0.1;
    Mesh out;
    // Reuse input vertices and intern new ones by position (bit-exact key).
    struct Vec3Key {
        double x, y, z;
        bool operator==(const Vec3Key& o) const {
            return x == o.x && y == o.y && z == o.z;
        }
    };
    struct Vec3Hash {
        size_t operator()(const Vec3Key& v) const {
            auto h1 = std::hash<double>{}(v.x);
            auto h2 = std::hash<double>{}(v.y);
            auto h3 = std::hash<double>{}(v.z);
            return h1 ^ (h2 * 16777619u) ^ (h3 * 2654435761u);
        }
    };
    std::unordered_map<Vec3Key, uint32_t, Vec3Hash> vmap;
    vmap.reserve(mesh.vertices.size() * 2);
    out.vertices.reserve(mesh.vertices.size() + result.pairs_found * 2);
    for (const auto& v : mesh.vertices) {
        Vec3Key k{ v[0], v[1], v[2] };
        vmap.emplace(k, static_cast<uint32_t>(out.vertices.size()));
        out.vertices.push_back(v);
    }
    auto intern = [&](const Vec3d& p) -> uint32_t {
        Vec3Key k{ p[0], p[1], p[2] };
        auto it = vmap.find(k);
        if (it != vmap.end()) return it->second;
        const uint32_t id = static_cast<uint32_t>(out.vertices.size());
        vmap.emplace(k, id);
        out.vertices.push_back({ p[0], p[1], p[2] });
        return id;
    };

    for (size_t fi = 0; fi < nf; ++fi) {
        const auto& f = F[fi];
        const Vec3d V0 = V[f[0]], V1 = V[f[1]], V2 = V[f[2]];
        if (tri_cuts[fi].empty()) {
            out.faces.push_back({ f[0], f[1], f[2] });
            continue;
        }
        auto cuts = consolidate(tri_cuts[fi], line_tol);
        auto pieces = internal::split_triangle_by_cuts(V0, V1, V2, cuts, edge_tol);
        if (pieces.empty()) {
            // Split produced nothing usable — keep original face.
            out.faces.push_back({ f[0], f[1], f[2] });
            continue;
        }
        for (const auto& p : pieces) {
            // Reject zero-area slivers from split (would be dropped by
            // degenerate stage anyway, but avoid emitting them).
            const Vec3d e1 = internal::sub(p[1], p[0]);
            const Vec3d e2 = internal::sub(p[2], p[0]);
            const Vec3d n  = internal::cross(e1, e2);
            if (internal::norm(n) < 1e-20) continue;
            const uint32_t a = intern(p[0]);
            const uint32_t b = intern(p[1]);
            const uint32_t c = intern(p[2]);
            if (a == b || b == c || a == c) continue;
            out.faces.push_back({ a, b, c });
        }
    }

    result.mesh         = std::move(out);
    result.faces_after  = static_cast<uint32_t>(result.mesh.faces.size());
    return result;
}

} // namespace meshseal::stages
