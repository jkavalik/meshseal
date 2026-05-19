#include "bridge_loops.h"
#include "../internal/diagnostics.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <unordered_map>
#include <vector>

namespace meshseal::stages {

namespace {

using Vec3 = std::array<double, 3>;

Vec3 sub(const Vec3& a, const Vec3& b) { return {a[0]-b[0], a[1]-b[1], a[2]-b[2]}; }
Vec3 cross(const Vec3& a, const Vec3& b) {
    return {a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]};
}
double dot(const Vec3& a, const Vec3& b) { return a[0]*b[0]+a[1]*b[1]+a[2]*b[2]; }
double dist(const Vec3& a, const Vec3& b) {
    const double dx=a[0]-b[0], dy=a[1]-b[1], dz=a[2]-b[2];
    return std::sqrt(dx*dx+dy*dy+dz*dz);
}

// --- union-find over vertices, grouping faces into edge-connected components.
struct DSU {
    std::vector<uint32_t> parent;
    explicit DSU(size_t n) : parent(n) {
        for (size_t i = 0; i < n; ++i) parent[i] = static_cast<uint32_t>(i);
    }
    uint32_t find(uint32_t x) {
        while (parent[x] != x) { parent[x] = parent[parent[x]]; x = parent[x]; }
        return x;
    }
    void unite(uint32_t a, uint32_t b) {
        uint32_t ra = find(a), rb = find(b);
        if (ra != rb) parent[ra] = rb;
    }
};

uint64_t ekey(uint32_t a, uint32_t b) {
    if (a > b) std::swap(a, b);
    return (static_cast<uint64_t>(a) << 32) | b;
}

// Trace open boundary loops (boundary edge = used by exactly one face).
// Returns vertex-index sequences in traversal order.
std::vector<std::vector<uint32_t>> boundary_loops(const Mesh& mesh) {
    std::unordered_map<uint64_t, int> count;
    std::unordered_map<uint64_t, std::array<uint32_t,2>> directed;
    count.reserve(mesh.faces.size() * 3);
    directed.reserve(mesh.faces.size() * 3);
    for (const auto& f : mesh.faces) {
        for (int i = 0; i < 3; ++i) {
            uint32_t a = f[i], b = f[(i+1)%3];
            uint64_t k = ekey(a, b);
            count[k]++;
            directed[k] = {a, b};
        }
    }
    std::unordered_map<uint32_t, std::vector<uint32_t>> adj;
    for (const auto& kv : count) {
        if (kv.second != 1) continue;
        const auto& d = directed[kv.first];
        adj[d[0]].push_back(d[1]);
        adj[d[1]].push_back(d[0]);
    }
    for (auto& kv : adj) std::sort(kv.second.begin(), kv.second.end());

    std::vector<std::vector<uint32_t>> loops;
    std::unordered_map<uint32_t, bool> visited;
    const uint32_t NONE = std::numeric_limits<uint32_t>::max();
    for (const auto& kv : adj) {
        uint32_t start = kv.first;
        if (visited[start]) continue;
        std::vector<uint32_t> loop{start};
        visited[start] = true;
        uint32_t prev = NONE, cur = start;
        for (;;) {
            uint32_t next = NONE;
            for (uint32_t n : adj[cur]) {
                if (!visited[n] && n != prev) { next = n; break; }
            }
            if (next == NONE) break;
            visited[next] = true;
            loop.push_back(next);
            prev = cur; cur = next;
        }
        if (loop.size() >= 3) loops.push_back(std::move(loop));
    }
    return loops;
}

// Möller–Trumbore ray-triangle intersection. Returns true and sets t (in units
// of |dir|) when the ray from orig along dir hits the triangle at t > eps.
bool ray_triangle(const Vec3& orig, const Vec3& dir,
                  const Vec3& v0, const Vec3& v1, const Vec3& v2, double& t) {
    const double eps = 1e-9;
    Vec3 e1 = sub(v1, v0), e2 = sub(v2, v0);
    Vec3 h = cross(dir, e2);
    double a = dot(e1, h);
    if (a > -eps && a < eps) return false;
    double f = 1.0 / a;
    Vec3 s = sub(orig, v0);
    double u = f * dot(s, h);
    if (u < 0.0 || u > 1.0) return false;
    Vec3 q = cross(s, e1);
    double v = f * dot(dir, q);
    if (v < 0.0 || u + v > 1.0) return false;
    t = f * dot(e2, q);
    return t > eps;
}

// True when the open segment a→b (with a tolerance band at both ends) is
// pierced by any face — i.e. solid geometry lies between the two centroids.
bool segment_blocked(const Mesh& mesh, const Vec3& a, const Vec3& b) {
    Vec3 d = sub(b, a);
    double seg_len = std::sqrt(dot(d, d));
    if (seg_len < 1e-12) return false;
    Vec3 unit{d[0]/seg_len, d[1]/seg_len, d[2]/seg_len};
    double tol = std::max(1.0, seg_len * 0.01);
    for (const auto& f : mesh.faces) {
        double t;
        if (ray_triangle(a, unit, mesh.vertices[f[0]], mesh.vertices[f[1]],
                         mesh.vertices[f[2]], t)) {
            if (t > tol && t < seg_len - tol) return true;
        }
    }
    return false;
}

// Greedy zipper: build a triangle strip bridging two open boundary loops.
std::vector<std::array<uint32_t,3>> bridge_two_loops(
        const std::vector<Vec3>& verts,
        const std::vector<uint32_t>& loopA,
        std::vector<uint32_t> loopB) {
    const int ma = static_cast<int>(loopA.size());
    int mb = static_cast<int>(loopB.size());
    std::vector<std::array<uint32_t,3>> faces;
    if (ma == 0 || mb == 0) return faces;

    // Align loopB: rotate so the vertex closest to loopA[0] comes first.
    const Vec3& va0 = verts[loopA[0]];
    int best_shift = 0;
    double best_d2 = std::numeric_limits<double>::max();
    for (int k = 0; k < mb; ++k) {
        const Vec3& p = verts[loopB[k]];
        double dx=p[0]-va0[0], dy=p[1]-va0[1], dz=p[2]-va0[2];
        double d2 = dx*dx+dy*dy+dz*dz;
        if (d2 < best_d2) { best_d2 = d2; best_shift = k; }
    }
    {
        std::vector<uint32_t> aligned(mb);
        for (int k = 0; k < mb; ++k) aligned[k] = loopB[(k+best_shift)%mb];
        loopB = std::move(aligned);
    }

    // Possibly reverse loopB[1:] so it runs the same direction as loopA.
    if (mb > 2 && ma > 1) {
        double d_fwd = dist(verts[loopB[1%mb]], verts[loopA[1%ma]]);
        double d_rev = dist(verts[loopB[mb-1]], verts[loopA[1%ma]]);
        if (d_rev < d_fwd) {
            std::vector<uint32_t> rev(mb);
            rev[0] = loopB[0];
            for (int k = 1; k < mb; ++k) rev[k] = loopB[mb-k];
            loopB = std::move(rev);
        }
    }

    // Greedy advance: at each step take the arc giving the shorter new edge.
    int ia = 0, ib = 0;
    while (ia < ma || ib < mb) {
        uint32_t a0 = loopA[ia%ma];
        uint32_t b0 = loopB[ib%mb];
        uint32_t a1 = loopA[(ia+1)%ma];
        uint32_t b1 = loopB[(ib+1)%mb];
        if (ia >= ma) {
            faces.push_back({a0, b0, b1});
            ++ib;
        } else if (ib >= mb) {
            faces.push_back({a0, b0, a1});
            ++ia;
        } else {
            double cost_a = dist(verts[a1], verts[b0]);
            double cost_b = dist(verts[a0], verts[b1]);
            if (cost_a <= cost_b) { faces.push_back({a0, b0, a1}); ++ia; }
            else                  { faces.push_back({a0, b0, b1}); ++ib; }
        }
    }
    return faces;
}

} // anonymous namespace

BridgeLoopsResult bridge_paired_loops(const Mesh& mesh,
                                      double   max_dist_factor,
                                      double   min_radius_ratio,
                                      double   min_vert_ratio,
                                      uint32_t min_loop_verts,
                                      uint32_t max_loop_verts) {
    BridgeLoopsResult result;
    result.mesh = mesh;

    const size_t nv = mesh.vertices.size();
    if (nv == 0 || mesh.faces.empty()) return result;

    // --- edge-connected component label per vertex ---
    DSU dsu(nv);
    for (const auto& f : mesh.faces) {
        dsu.unite(f[0], f[1]);
        dsu.unite(f[1], f[2]);
    }

    // --- candidate boundary loops ---
    auto loops = boundary_loops(mesh);
    struct LoopInfo {
        std::vector<uint32_t> loop;
        Vec3     centroid;
        double   radius;
        uint32_t comp;
    };
    std::vector<LoopInfo> infos;
    for (auto& loop : loops) {
        if (loop.size() < min_loop_verts || loop.size() > max_loop_verts) continue;
        Vec3 c{0,0,0};
        for (uint32_t vi : loop) {
            c[0] += mesh.vertices[vi][0];
            c[1] += mesh.vertices[vi][1];
            c[2] += mesh.vertices[vi][2];
        }
        double n = static_cast<double>(loop.size());
        c[0] /= n; c[1] /= n; c[2] /= n;
        double max_r = 0.0;
        for (uint32_t vi : loop) {
            double r = dist(mesh.vertices[vi], c);
            if (r > max_r) max_r = r;
        }
        uint32_t comp = dsu.find(loop[0]);
        infos.push_back({std::move(loop), c, max_r, comp});
    }
    if (infos.size() < 2) return result;

    // --- cross-component pairs sorted by centroid distance ---
    struct PairEntry { double d; int i, j; };
    std::vector<PairEntry> pairs;
    for (size_t i = 0; i < infos.size(); ++i)
        for (size_t j = i+1; j < infos.size(); ++j) {
            if (infos[i].comp == infos[j].comp) continue;
            double d = dist(infos[i].centroid, infos[j].centroid);
            pairs.push_back({d, static_cast<int>(i), static_cast<int>(j)});
        }
    std::sort(pairs.begin(), pairs.end(),
              [](const PairEntry& a, const PairEntry& b){ return a.d < b.d; });

    // --- greedy bridge selection: each loop participates at most once ---
    // Each accepted strip is guarded: a bridge that does not strictly REDUCE
    // the open-boundary count, or that INCREASES the non-manifold-edge count,
    // is a mis-pair (the two loops were not a genuine shared opening) — its
    // strip is rolled back. This keeps captain_toad's mouth bridges while
    // rejecting the spurious cross-component pairings on bird_bath / human /
    // sea_vase / black_vase.
    std::vector<bool> bridged(infos.size(), false);
    std::vector<Vec3> verts(mesh.vertices.begin(), mesh.vertices.end());
    auto base_diag = internal::compute_diagnostics(result.mesh);
    int  cur_nm  = static_cast<int>(base_diag.non_manifold_edges);
    int  cur_bnd = static_cast<int>(base_diag.open_boundary_edges);
    for (const auto& pe : pairs) {
        int i = pe.i, j = pe.j;
        if (bridged[i] || bridged[j]) continue;
        double avg_r = (infos[i].radius + infos[j].radius) * 0.5;
        if (avg_r < 1e-9 || pe.d >= max_dist_factor * avg_r) continue;
        // Genuinely-coincident gate: similar radius and similar vertex count.
        double ri = infos[i].radius, rj = infos[j].radius;
        double r_ratio = std::min(ri, rj) / std::max(ri, rj);
        if (r_ratio < min_radius_ratio) continue;
        double ni = static_cast<double>(infos[i].loop.size());
        double nj = static_cast<double>(infos[j].loop.size());
        double v_ratio = std::min(ni, nj) / std::max(ni, nj);
        if (v_ratio < min_vert_ratio) continue;
        if (segment_blocked(mesh, infos[i].centroid, infos[j].centroid)) continue;

        auto bf = bridge_two_loops(verts, infos[i].loop, infos[j].loop);
        if (bf.empty()) continue;

        const size_t before = result.mesh.faces.size();
        for (const auto& t : bf) result.mesh.faces.push_back(t);
        auto post = internal::compute_diagnostics(result.mesh);
        const int post_nm  = static_cast<int>(post.non_manifold_edges);
        const int post_bnd = static_cast<int>(post.open_boundary_edges);
        if (post_nm > cur_nm || post_bnd >= cur_bnd) {
            // Mis-pair — roll the strip back.
            result.mesh.faces.resize(before);
            continue;
        }
        cur_nm = post_nm;
        cur_bnd = post_bnd;
        result.faces_added += static_cast<uint32_t>(bf.size());
        bridged[i] = bridged[j] = true;
        ++result.bridges_made;
    }
    return result;
}

} // namespace meshseal::stages
