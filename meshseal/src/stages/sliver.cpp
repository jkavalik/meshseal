#include "sliver.h"
#include "../internal/vec3.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <unordered_map>
#include <vector>

namespace meshseal::stages {

using internal::Vec3d;

namespace {

uint64_t ekey(uint32_t a, uint32_t b) {
    if (a > b) std::swap(a, b);
    return (static_cast<uint64_t>(a) << 32) | b;
}

// Unnormalised face normal (cross product) — magnitude is 2·area.
Vec3d face_cross(const Vec3d& a, const Vec3d& b, const Vec3d& c) {
    return internal::cross(internal::sub(b, a), internal::sub(c, a));
}

double tri_quality(const Vec3d& a, const Vec3d& b, const Vec3d& c) {
    const Vec3d n = face_cross(a, b, c);
    const double area2 = internal::norm(n);            // = 2·area
    const Vec3d e0 = internal::sub(b, a);
    const Vec3d e1 = internal::sub(c, b);
    const Vec3d e2 = internal::sub(a, c);
    const double sum_sq = internal::dot(e0, e0) +
                          internal::dot(e1, e1) +
                          internal::dot(e2, e2);
    if (sum_sq <= 0.0) return 0.0;
    // q = 4√3·area / Σℓ²  = 2√3·(2·area) / Σℓ²
    return (2.0 * 1.7320508075688772 * area2) / sum_sq;
}

} // anonymous namespace

SliverResult collapse_slivers(const Mesh& mesh, double quality_threshold) {
    SliverResult result;
    result.mesh = mesh;

    const uint32_t nverts = static_cast<uint32_t>(mesh.vertices.size());
    if (mesh.faces.empty() || nverts == 0) return result;

    // Working copy.
    std::vector<std::array<uint32_t, 3>> faces(mesh.faces.begin(), mesh.faces.end());
    std::vector<Vec3d> verts = mesh.vertices;
    std::vector<bool> face_dead(faces.size(), false);

    // Passes: each pass collapses a set of non-conflicting slivers (no two
    // collapses sharing a vertex), then repeats. Bounded — a collapse can
    // create a new sliver, but each collapse strictly removes ≥2 faces, so
    // the face count monotonically decreases; cap at a generous pass limit.
    const int kMaxPasses = 12;
    for (int pass = 0; pass < kMaxPasses; ++pass) {
        // Build edge → incident (live) faces.
        std::unordered_map<uint64_t, std::vector<uint32_t>> e2f;
        e2f.reserve(faces.size() * 3);
        for (uint32_t fi = 0; fi < faces.size(); ++fi) {
            if (face_dead[fi]) continue;
            const auto& f = faces[fi];
            e2f[ekey(f[0], f[1])].push_back(fi);
            e2f[ekey(f[1], f[2])].push_back(fi);
            e2f[ekey(f[2], f[0])].push_back(fi);
        }
        // vertex → incident (live) faces.
        std::vector<std::vector<uint32_t>> v2f(nverts);
        for (uint32_t fi = 0; fi < faces.size(); ++fi) {
            if (face_dead[fi]) continue;
            for (int k = 0; k < 3; ++k) v2f[faces[fi][k]].push_back(fi);
        }

        // Collect sliver faces, worst quality first.
        std::vector<std::pair<double, uint32_t>> slivers;
        for (uint32_t fi = 0; fi < faces.size(); ++fi) {
            if (face_dead[fi]) continue;
            const auto& f = faces[fi];
            const double q = tri_quality(verts[f[0]], verts[f[1]], verts[f[2]]);
            if (q < quality_threshold) slivers.emplace_back(q, fi);
        }
        if (slivers.empty()) break;
        std::sort(slivers.begin(), slivers.end());

        std::vector<bool> vtouched(nverts, false);
        bool any = false;

        for (const auto& [q, fi] : slivers) {
            (void)q;
            if (face_dead[fi]) continue;
            const auto& f = faces[fi];

            // Shortest edge of this sliver.
            int best_e = -1;
            double best_len = 0.0;
            for (int k = 0; k < 3; ++k) {
                const Vec3d d = internal::sub(verts[f[(k+1)%3]], verts[f[k]]);
                const double l = internal::norm(d);
                if (best_e < 0 || l < best_len) { best_len = l; best_e = k; }
            }
            uint32_t u = f[best_e];
            uint32_t v = f[(best_e + 1) % 3];
            if (u == v) continue;
            if (vtouched[u] || vtouched[v]) continue;  // conflict — defer to next pass

            // The shortest edge must be MANIFOLD (exactly 2 faces). Collapsing
            // a non-manifold or boundary short edge tears the surface.
            auto it = e2f.find(ekey(u, v));
            if (it == e2f.end() || it->second.size() != 2) continue;

            // Tentative collapse v → u. Validate guards.
            bool ok = true;
            // Faces incident to v that survive (don't also contain u):
            for (uint32_t nf : v2f[v]) {
                if (face_dead[nf]) continue;
                const auto& g = faces[nf];
                bool has_u = (g[0]==u || g[1]==u || g[2]==u);
                if (has_u) continue;  // this face shares edge (u,v) → will die
                // Build post-collapse face (v→u).
                std::array<uint32_t,3> ng = g;
                for (int k = 0; k < 3; ++k) if (ng[k]==v) ng[k]=u;
                if (ng[0]==ng[1] || ng[1]==ng[2] || ng[2]==ng[0]) {
                    // Would become degenerate even though it didn't contain u
                    // pre-collapse — impossible unless mesh is malformed; reject.
                    ok = false; break;
                }
                const Vec3d n_pre  = face_cross(verts[g[0]],  verts[g[1]],  verts[g[2]]);
                const Vec3d n_post = face_cross(verts[ng[0]], verts[ng[1]], verts[ng[2]]);
                if (internal::dot(n_pre, n_post) <= 0.0) { ok = false; break; }  // normal flip
            }
            if (!ok) continue;

            // Non-manifold creation guard: an edge (v,x) collapsing onto an
            // existing (u,x) must not push the merged edge past 2 faces.
            for (uint32_t nf : v2f[v]) {
                if (face_dead[nf]) continue;
                const auto& g = faces[nf];
                for (int k = 0; k < 3 && ok; ++k) {
                    uint32_t a = g[k], b = g[(k+1)%3];
                    if (a != v && b != v) continue;
                    uint32_t x = (a == v) ? b : a;
                    if (x == u) continue;  // the collapsing edge itself
                    // post-collapse this edge becomes (u,x). Count live faces
                    // currently on (u,x) and on (v,x); the (u,v)-pair faces die.
                    auto iu = e2f.find(ekey(u, x));
                    auto iv = e2f.find(ekey(v, x));
                    int cu = (iu==e2f.end()) ? 0 : static_cast<int>(iu->second.size());
                    int cv = (iv==e2f.end()) ? 0 : static_cast<int>(iv->second.size());
                    if (cu + cv > 2) { ok = false; }
                }
            }
            if (!ok) continue;

            // Apply: v → u everywhere; kill faces that become degenerate.
            for (uint32_t nf : v2f[v]) {
                if (face_dead[nf]) continue;
                auto& g = faces[nf];
                for (int k = 0; k < 3; ++k) if (g[k]==v) g[k]=u;
                if (g[0]==g[1] || g[1]==g[2] || g[2]==g[0]) face_dead[nf] = true;
            }
            vtouched[u] = true;
            vtouched[v] = true;
            ++result.slivers_collapsed;
            any = true;
        }
        if (!any) break;
    }

    // Rebuild output mesh: live faces, compacted vertices.
    Mesh out;
    std::vector<uint32_t> remap(nverts, UINT32_MAX);
    for (uint32_t fi = 0; fi < faces.size(); ++fi) {
        if (face_dead[fi]) continue;
        std::array<uint32_t,3> nf{};
        for (int k = 0; k < 3; ++k) {
            uint32_t ov = faces[fi][k];
            if (remap[ov] == UINT32_MAX) {
                remap[ov] = static_cast<uint32_t>(out.vertices.size());
                out.vertices.push_back(verts[ov]);
            }
            nf[k] = remap[ov];
        }
        out.faces.push_back(nf);
    }
    result.mesh = std::move(out);
    return result;
}

} // namespace meshseal::stages
