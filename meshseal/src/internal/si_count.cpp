#include "si_count.h"
#include "arrangement.h"
#include "vec3.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace meshseal::internal {

namespace {

static uint64_t cell_key(int32_t x, int32_t y, int32_t z) {
    auto pack = [](int32_t v) -> uint64_t {
        return static_cast<uint64_t>(static_cast<int64_t>(v) + (1 << 20)) & 0x1FFFFF;
    };
    return pack(x) | (pack(y) << 21) | (pack(z) << 42);
}

static uint64_t ekey(uint32_t a, uint32_t b) {
    if (a > b) std::swap(a, b);
    return (static_cast<uint64_t>(a) << 32) | b;
}

} // anonymous namespace

int count_self_intersections(const Mesh& mesh) {
    const std::size_t nf = mesh.faces.size();
    if (nf < 2 || mesh.vertices.empty()) return 0;
    // Skip on very large meshes — even with bucket pruning, O(F * candidates)
    // dominates total repair time. Treat as "not measured" rather than
    // producing a misleading 0 or burning minutes on dense scan inputs.
    if (nf > 40000) return -1;
    // Cap the accumulated count at 10000 — the user already knows there's
    // a problem once we hit this many pairs; further counting wastes time.
    constexpr int kCountCap = 10000;

    using Vec3 = std::array<double, 3>;
    const auto& V = mesh.vertices;
    const auto& F = mesh.faces;

    // bbox + median edge length → cell size
    Vec3 lo = V[0], hi = V[0];
    for (const auto& p : V)
        for (int k = 0; k < 3; ++k) {
            lo[k] = std::min(lo[k], p[k]);
            hi[k] = std::max(hi[k], p[k]);
        }
    const double dx = hi[0]-lo[0], dy = hi[1]-lo[1], dz = hi[2]-lo[2];
    const double bbox_diag = std::sqrt(dx*dx + dy*dy + dz*dz);
    if (bbox_diag < 1e-30) return 0;

    std::vector<double> edge_lens;
    edge_lens.reserve(nf * 3);
    for (const auto& t : F) {
        for (int k = 0; k < 3; ++k) {
            const auto& a = V[t[k]]; const auto& b = V[t[(k+1)%3]];
            const double ex = a[0]-b[0], ey = a[1]-b[1], ez = a[2]-b[2];
            edge_lens.push_back(std::sqrt(ex*ex + ey*ey + ez*ez));
        }
    }
    std::nth_element(edge_lens.begin(),
                     edge_lens.begin() + edge_lens.size()/2,
                     edge_lens.end());
    double cell = edge_lens[edge_lens.size()/2];
    if (cell < bbox_diag * 1e-6) cell = bbox_diag * 1e-6;

    // Bucket triangles by their bbox cells.
    std::unordered_map<uint64_t, std::vector<uint32_t>> grid;
    grid.reserve(nf * 2);
    std::vector<std::array<Vec3,2>> tri_bb(nf);
    for (std::size_t fi = 0; fi < nf; ++fi) {
        const auto& t = F[fi];
        Vec3 mn{ V[t[0]][0], V[t[0]][1], V[t[0]][2] }, mx = mn;
        for (int k = 1; k < 3; ++k) {
            for (int j = 0; j < 3; ++j) {
                mn[j] = std::min(mn[j], V[t[k]][j]);
                mx[j] = std::max(mx[j], V[t[k]][j]);
            }
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

    // Edge-adjacency: skip pairs that share an edge (manifold adjacency
    // produces a degenerate "intersection" along the shared edge).
    std::unordered_map<uint64_t, std::vector<uint32_t>> e2f;
    e2f.reserve(nf * 3);
    for (std::size_t fi = 0; fi < nf; ++fi) {
        const auto& t = F[fi];
        for (int k = 0; k < 3; ++k)
            e2f[ekey(t[k], t[(k+1)%3])].push_back(static_cast<uint32_t>(fi));
    }
    // Skip any pair sharing one or more vertices. Fan-around-vertex pairs
    // (1 shared) and edge-adjacent pairs (2 shared) can give Möller a
    // numerical-artifact "intersection segment" at the shared vertex/edge
    // that is not a genuine self-intersection. Matches the convention of
    // the external `si_count.py` diagnostic.
    auto shares_vertex = [&](uint32_t i, uint32_t j) {
        const auto& ti = F[i]; const auto& tj = F[j];
        for (int a = 0; a < 3; ++a)
            for (int b = 0; b < 3; ++b)
                if (ti[a] == tj[b]) return true;
        return false;
    };

    auto bb_overlap = [](const std::array<Vec3,2>& a, const std::array<Vec3,2>& b) {
        for (int k = 0; k < 3; ++k)
            if (a[1][k] < b[0][k] || b[1][k] < a[0][k]) return false;
        return true;
    };

    // For each face, collect candidate partners from the grid; T-T test each.
    std::unordered_set<uint64_t> hits;
    hits.reserve(nf);
    for (std::size_t fi = 0; fi < nf; ++fi) {
        const auto& mn = tri_bb[fi][0]; const auto& mx = tri_bb[fi][1];
        const int32_t cx0 = static_cast<int32_t>(std::floor(mn[0] / cell));
        const int32_t cy0 = static_cast<int32_t>(std::floor(mn[1] / cell));
        const int32_t cz0 = static_cast<int32_t>(std::floor(mn[2] / cell));
        const int32_t cx1 = static_cast<int32_t>(std::floor(mx[0] / cell));
        const int32_t cy1 = static_cast<int32_t>(std::floor(mx[1] / cell));
        const int32_t cz1 = static_cast<int32_t>(std::floor(mx[2] / cell));
        std::unordered_set<uint32_t> cands;
        for (int32_t cx = cx0; cx <= cx1; ++cx)
        for (int32_t cy = cy0; cy <= cy1; ++cy)
        for (int32_t cz = cz0; cz <= cz1; ++cz) {
            auto it = grid.find(cell_key(cx, cy, cz));
            if (it == grid.end()) continue;
            for (uint32_t fj : it->second) if (fj > fi) cands.insert(fj);
        }
        for (uint32_t fj : cands) {
            if (!bb_overlap(tri_bb[fi], tri_bb[fj])) continue;
            if (shares_vertex(static_cast<uint32_t>(fi), fj)) continue;
            const auto& ti = F[fi]; const auto& tj = F[fj];
            // Convert to Vec3d for tri_tri_intersect.
            const Vec3d a1 = V[ti[0]], b1 = V[ti[1]], c1 = V[ti[2]];
            const Vec3d a2 = V[tj[0]], b2 = V[tj[1]], c2 = V[tj[2]];
            auto seg = tri_tri_intersect(a1, b1, c1, a2, b2, c2);
            if (seg.size() < 2) continue;   // empty / coplanar / degenerate
            // Ignore shared-vertex "intersections" that produced a zero-
            // length or tiny segment at the vertex itself.
            const double ddx = seg[1][0] - seg[0][0];
            const double ddy = seg[1][1] - seg[0][1];
            const double ddz = seg[1][2] - seg[0][2];
            const double len2 = ddx*ddx + ddy*ddy + ddz*ddz;
            if (len2 < bbox_diag * bbox_diag * 1e-20) continue;
            uint64_t key = (static_cast<uint64_t>(fi) << 32) | fj;
            hits.insert(key);
            if (static_cast<int>(hits.size()) >= kCountCap) goto done;
        }
    }
done:
    return static_cast<int>(hits.size());
}

} // namespace meshseal::internal
