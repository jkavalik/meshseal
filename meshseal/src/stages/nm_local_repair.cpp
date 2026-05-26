#include "nm_local_repair.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace meshseal::stages {

namespace {

using Tri  = std::array<uint32_t, 3>;
using Vec3 = std::array<double, 3>;

static uint64_t ekey(uint32_t a, uint32_t b) {
    if (a > b) std::swap(a, b);
    return (static_cast<uint64_t>(a) << 32) | b;
}

static double dist2(const Vec3& a, const Vec3& b) {
    const double dx = a[0]-b[0], dy = a[1]-b[1], dz = a[2]-b[2];
    return dx*dx + dy*dy + dz*dz;
}

// Canonical key for a face's vertex set: sorted triple (a,b,c) with a<b<c.
struct VSet {
    uint32_t a, b, c;
    bool operator==(const VSet& o) const { return a==o.a && b==o.b && c==o.c; }
};
struct VSetHash {
    size_t operator()(const VSet& v) const {
        return (size_t(v.a)*73856093u) ^ (size_t(v.b)*19349663u) ^ (size_t(v.c)*83492791u);
    }
};
static VSet vset(const Tri& t) {
    uint32_t a = t[0], b = t[1], c = t[2];
    if (a > b) std::swap(a, b);
    if (b > c) std::swap(b, c);
    if (a > b) std::swap(a, b);
    return {a, b, c};
}

// Winding sign of a face relative to its sorted-canonical form. Two
// faces with the same vset have OPPOSITE winding iff their winding
// signs differ.
static int winding_sign(const Tri& t) {
    int m = 0;
    if (t[1] < t[m]) m = 1;
    if (t[2] < t[m]) m = 2;
    uint32_t a = t[m];
    uint32_t b = t[(m+1) % 3];
    uint32_t c = t[(m+2) % 3];
    (void)a;
    return (b < c) ? 1 : -1;
}

} // anonymous namespace

NmLocalRepairResult nm_local_repair(const Mesh& mesh, double tol_rel) {
    NmLocalRepairResult result;
    result.mesh = mesh;

    const uint32_t nf = static_cast<uint32_t>(mesh.faces.size());
    if (nf == 0 || mesh.vertices.empty()) return result;

    // bbox diag → absolute weld tolerance
    Vec3 lo = mesh.vertices[0], hi = mesh.vertices[0];
    for (const auto& p : mesh.vertices)
        for (int k = 0; k < 3; ++k) {
            lo[k] = std::min(lo[k], p[k]);
            hi[k] = std::max(hi[k], p[k]);
        }
    const double dx = hi[0]-lo[0], dy = hi[1]-lo[1], dz = hi[2]-lo[2];
    const double bbox_diag = std::sqrt(dx*dx + dy*dy + dz*dz);
    const double tol = bbox_diag * tol_rel;
    const double tol2 = tol * tol;

    // ---------- Phase 1: NM-local proximity weld ----------
    // edge -> incident face list; vertex -> 1-ring vertex set
    std::unordered_map<uint64_t, std::vector<uint32_t>> e2f;
    e2f.reserve(nf * 3);
    std::vector<std::unordered_set<uint32_t>> v_neigh(mesh.vertices.size());
    for (uint32_t fi = 0; fi < nf; ++fi) {
        const auto& f = mesh.faces[fi];
        for (int k = 0; k < 3; ++k) {
            const uint32_t a = f[k];
            const uint32_t b = f[(k+1) % 3];
            e2f[ekey(a, b)].push_back(fi);
            v_neigh[a].insert(b);
            v_neigh[b].insert(a);
        }
    }

    // Collect NM-edge endpoints
    std::unordered_set<uint32_t> nm_endpoints;
    for (const auto& kv : e2f) {
        if (kv.second.size() > 2) {
            uint32_t a = static_cast<uint32_t>(kv.first >> 32);
            uint32_t b = static_cast<uint32_t>(kv.first & 0xffffffffu);
            nm_endpoints.insert(a);
            nm_endpoints.insert(b);
        }
    }
    if (nm_endpoints.empty()) return result;

    // NOTE: iteration is unordered_set-bucket-defined. Documented as a
    // non-determinism risk. Sorting REGRESSES black_vase (rescue depends
    // on specific merge survivor choice). Left as-is until downstream
    // decoupling. See REVIEW_CHECKLIST Phase 3 reversion note.

    // Union-find for vertex merges. Always merge w → ep (the NM endpoint
    // is the survivor; w is the duplicate that gets absorbed).
    std::vector<uint32_t> uf(mesh.vertices.size());
    for (uint32_t i = 0; i < uf.size(); ++i) uf[i] = i;
    auto find = [&](uint32_t x) {
        while (uf[x] != x) { uf[x] = uf[uf[x]]; x = uf[x]; }
        return x;
    };
    uint32_t merges = 0;
    for (uint32_t ep : nm_endpoints) {
        for (uint32_t w : v_neigh[ep]) {
            if (w == ep) continue;
            const uint32_t r_ep = find(ep);
            const uint32_t r_w  = find(w);
            if (r_ep == r_w) continue;
            if (dist2(mesh.vertices[r_ep], mesh.vertices[r_w]) < tol2) {
                uf[r_w] = r_ep;
                ++merges;
            }
        }
    }
    result.merges = merges;

    // Apply union-find to faces, drop degenerates
    std::vector<Tri> faces1;
    faces1.reserve(nf);
    for (const auto& f : mesh.faces) {
        const uint32_t a = find(f[0]), b = find(f[1]), c = find(f[2]);
        if (a == b || b == c || a == c) continue;
        faces1.push_back({a, b, c});
    }

    // ---------- Phase 2: strict back-to-back pair dedup ----------
    std::unordered_map<VSet, std::vector<uint32_t>, VSetHash> by_set;
    by_set.reserve(faces1.size());
    for (uint32_t fi = 0; fi < faces1.size(); ++fi) by_set[vset(faces1[fi])].push_back(fi);
    std::vector<bool> remove(faces1.size(), false);
    uint32_t pairs = 0;
    for (auto& kv : by_set) {
        const auto& lst = kv.second;
        if (lst.size() < 2) continue;
        // Greedy: pair the first +winding with the first -winding, etc.
        std::vector<uint32_t> plus, minus;
        for (uint32_t fi : lst) {
            if (remove[fi]) continue;
            (winding_sign(faces1[fi]) > 0 ? plus : minus).push_back(fi);
        }
        const size_t n = std::min(plus.size(), minus.size());
        for (size_t i = 0; i < n; ++i) {
            remove[plus[i]]  = true;
            remove[minus[i]] = true;
            ++pairs;
        }
    }
    result.pairs_removed = pairs;

    std::vector<Tri> faces2;
    faces2.reserve(faces1.size());
    for (uint32_t fi = 0; fi < faces1.size(); ++fi)
        if (!remove[fi]) faces2.push_back(faces1[fi]);

    result.faces_dropped = nf - static_cast<uint32_t>(faces2.size());

    // Build output mesh — keep vertices array, faces only via remap of
    // the union-find representatives (unreferenced vertices harmless).
    Mesh out;
    out.vertices = mesh.vertices;
    out.faces.reserve(faces2.size());
    for (const auto& t : faces2) out.faces.push_back({t[0], t[1], t[2]});
    result.mesh = std::move(out);
    return result;
}

} // namespace meshseal::stages
