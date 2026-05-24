#include "coplanar_fan_drop.h"
#include <array>
#include <cmath>
#include <cstdint>
#include <unordered_map>
#include <vector>

namespace meshseal::stages {

namespace {

inline uint64_t ekey(uint32_t a, uint32_t b) {
    if (a > b) std::swap(a, b);
    return (static_cast<uint64_t>(a) << 32) | b;
}

struct Vec3d { double x, y, z; };

inline Vec3d cross(const Vec3d& a, const Vec3d& b) {
    return { a.y * b.z - a.z * b.y,
             a.z * b.x - a.x * b.z,
             a.x * b.y - a.y * b.x };
}

inline double dot(const Vec3d& a, const Vec3d& b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

inline double norm(const Vec3d& a) {
    return std::sqrt(dot(a, a));
}

inline Vec3d sub(const std::array<double, 3>& p, const std::array<double, 3>& q) {
    return { p[0] - q[0], p[1] - q[1], p[2] - q[2] };
}

struct FaceGeo {
    Vec3d  normal;   // unit, or {0,0,0} if degenerate
    double area;
};

inline FaceGeo face_geom(const Mesh& m, uint32_t fi) {
    const auto& f = m.faces[fi];
    Vec3d e1 = sub(m.vertices[f[1]], m.vertices[f[0]]);
    Vec3d e2 = sub(m.vertices[f[2]], m.vertices[f[0]]);
    Vec3d n  = cross(e1, e2);
    double L = norm(n);
    if (L < 1e-30) return { {0,0,0}, 0.0 };
    return { { n.x / L, n.y / L, n.z / L }, L * 0.5 };
}

} // namespace

CoplanarFanDropResult coplanar_fan_drop(const Mesh& mesh, double coplanar_dot) {
    CoplanarFanDropResult result;
    result.mesh = mesh;
    const uint32_t nf = static_cast<uint32_t>(mesh.faces.size());
    if (nf < 4) return result;

    // Edge -> list of (face_index, directed_a, directed_b) so we can
    // distinguish forward vs reverse winding on the edge.
    struct Inc { uint32_t fi; uint32_t a; uint32_t b; };
    std::unordered_map<uint64_t, std::vector<Inc>> e2f;
    e2f.reserve(nf * 2);
    for (uint32_t fi = 0; fi < nf; ++fi) {
        const auto& f = mesh.faces[fi];
        for (int k = 0; k < 3; ++k) {
            uint32_t a = f[k], b = f[(k + 1) % 3];
            e2f[ekey(a, b)].push_back({fi, a, b});
        }
    }

    std::vector<bool> drop(nf, false);

    for (const auto& kv : e2f) {
        const auto& incs = kv.second;
        if (incs.size() != 4) continue;
        // canonical edge ordering (lo < hi)
        uint32_t lo = static_cast<uint32_t>(kv.first >> 32);
        uint32_t hi = static_cast<uint32_t>(kv.first & 0xffffffffu);
        // Group by winding direction
        uint32_t fwd[2] = {UINT32_MAX, UINT32_MAX};
        uint32_t rev[2] = {UINT32_MAX, UINT32_MAX};
        int nf_ = 0, nr_ = 0;
        for (const auto& inc : incs) {
            if (inc.a == lo && inc.b == hi) {
                if (nf_ < 2) fwd[nf_++] = inc.fi;
            } else {
                if (nr_ < 2) rev[nr_++] = inc.fi;
            }
        }
        if (nf_ != 2 || nr_ != 2) continue;  // not a 2+2 fan
        if (drop[fwd[0]] || drop[fwd[1]] || drop[rev[0]] || drop[rev[1]])
            continue;  // already touched by another edge's resolution

        // All 4 faces must be coplanar (within tolerance).
        FaceGeo g[4];
        g[0] = face_geom(mesh, fwd[0]);
        g[1] = face_geom(mesh, fwd[1]);
        g[2] = face_geom(mesh, rev[0]);
        g[3] = face_geom(mesh, rev[1]);
        bool all_coplanar = true;
        for (int i = 0; i < 4 && all_coplanar; ++i) {
            if (g[i].area <= 0) { all_coplanar = false; break; }
            for (int j = i + 1; j < 4; ++j) {
                if (std::abs(dot(g[i].normal, g[j].normal)) < coplanar_dot) {
                    all_coplanar = false; break;
                }
            }
        }
        if (!all_coplanar) continue;

        // Two manifold pairings: (fwd0-rev0, fwd1-rev1) or (fwd0-rev1, fwd1-rev0).
        // Drop the pair whose total area is smaller — that's the
        // "duplicate" layer.
        const double pa = g[0].area + g[2].area;   // fwd0+rev0
        const double pb = g[1].area + g[3].area;   // fwd1+rev1
        const double pc = g[0].area + g[3].area;   // fwd0+rev1
        const double pd = g[1].area + g[2].area;   // fwd1+rev0
        // Pairing 1: (a, b); pairing 2: (c, d). For each, the smaller
        // pair-area is the "drop candidate"; pick the pairing where the
        // ratio (smaller / larger) is most lopsided — most confident
        // identification of the duplicate layer.
        const double r1 = (std::min(pa, pb) / std::max(pa, pb));
        const double r2 = (std::min(pc, pd) / std::max(pc, pd));
        bool use_pairing1 = (r1 < r2);
        uint32_t drop1, drop2;
        if (use_pairing1) {
            if (pa < pb) { drop1 = fwd[0]; drop2 = rev[0]; }
            else         { drop1 = fwd[1]; drop2 = rev[1]; }
        } else {
            if (pc < pd) { drop1 = fwd[0]; drop2 = rev[1]; }
            else         { drop1 = fwd[1]; drop2 = rev[0]; }
        }
        drop[drop1] = true;
        drop[drop2] = true;
        ++result.fans_dropped;
        result.faces_removed += 2;
    }

    if (result.fans_dropped == 0) return result;

    // Build new mesh without dropped faces. Vertices kept as-is
    // (unreferenced vertices are harmless downstream).
    Mesh out;
    out.vertices = mesh.vertices;
    out.faces.reserve(nf - result.faces_removed);
    for (uint32_t fi = 0; fi < nf; ++fi) {
        if (!drop[fi]) out.faces.push_back(mesh.faces[fi]);
    }
    result.mesh    = std::move(out);
    result.applied = true;
    return result;
}

} // namespace meshseal::stages
