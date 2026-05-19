#include "collapse_nm.h"
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

uint64_t ekey(uint32_t a, uint32_t b) {
    if (a > b) std::swap(a, b);
    return (static_cast<uint64_t>(a) << 32) | b;
}
void edecode(uint64_t k, uint32_t& a, uint32_t& b) {
    a = static_cast<uint32_t>(k >> 32);
    b = static_cast<uint32_t>(k & 0xffffffffu);
}

double vdist(const Vec3& a, const Vec3& b) {
    double dx = a[0]-b[0], dy = a[1]-b[1], dz = a[2]-b[2];
    return std::sqrt(dx*dx + dy*dy + dz*dz);
}

// (B-A)x(C-A); magnitude == 2·area.
Vec3 tri_cross(const Vec3& A, const Vec3& B, const Vec3& C) {
    double ux = B[0]-A[0], uy = B[1]-A[1], uz = B[2]-A[2];
    double wx = C[0]-A[0], wy = C[1]-A[1], wz = C[2]-A[2];
    return { uy*wz - uz*wy, uz*wx - ux*wz, ux*wy - uy*wx };
}

// Signed volume of the tetrahedron (origin, A, B, C) times 6.
double tet6(const Vec3& A, const Vec3& B, const Vec3& C) {
    return A[0]*(B[1]*C[2]-B[2]*C[1])
         - A[1]*(B[0]*C[2]-B[2]*C[0])
         + A[2]*(B[0]*C[1]-B[1]*C[0]);
}

} // anonymous namespace

CollapseNmResult collapse_nm_region(const Mesh& mesh, double max_collapse_len,
                                    int max_collapses) {
    CollapseNmResult result;
    result.mesh = mesh;
    if (mesh.faces.empty() || mesh.vertices.empty()) return result;

    std::vector<Vec3> V(mesh.vertices.begin(), mesh.vertices.end());
    std::vector<Tri>  F(mesh.faces.begin(), mesh.faces.end());

    Vec3 lo = V[0], hi = V[0];
    for (const auto& p : V)
        for (int k = 0; k < 3; ++k) {
            lo[k] = std::min(lo[k], p[k]);
            hi[k] = std::max(hi[k], p[k]);
        }
    const double bbox_diag = vdist(lo, hi);
    // Only ever collapse SHORT edges — the slivers and pinch edges of a
    // degenerate membrane. Real feature edges (and the wall thickness of a
    // thin-walled model) stay safely above the cap.
    if (max_collapse_len < 0.0) max_collapse_len = bbox_diag * 0.006;

    // --- edge incidence map ---
    std::unordered_map<uint64_t, int> e2c;
    e2c.reserve(F.size() * 3);
    auto add_face = [&](const Tri& f) {
        for (int k = 0; k < 3; ++k) e2c[ekey(f[k], f[(k+1)%3])]++;
    };
    auto sub_face = [&](const Tri& f) {
        for (int k = 0; k < 3; ++k) {
            auto it = e2c.find(ekey(f[k], f[(k+1)%3]));
            if (it != e2c.end() && --it->second == 0) e2c.erase(it);
        }
    };
    for (const auto& f : F) add_face(f);

    int nm = 0, bnd = 0;
    for (const auto& kv : e2c) {
        if (kv.second > 2) ++nm;
        else if (kv.second == 1) ++bnd;
    }
    result.nm_before = static_cast<uint32_t>(nm);
    if (nm == 0) return result;
    result.applied = true;

    std::vector<char> alive(F.size(), 1);
    std::unordered_map<uint32_t, std::vector<uint32_t>> v2f;
    for (uint32_t fi = 0; fi < F.size(); ++fi)
        for (int k = 0; k < 3; ++k) v2f[F[fi][k]].push_back(fi);

    // sliver area threshold: 5 % of the median 2·area.
    double tau = 1e-12;
    {
        std::vector<double> mags;
        mags.reserve(F.size());
        for (const auto& f : F) {
            Vec3 c = tri_cross(V[f[0]], V[f[1]], V[f[2]]);
            mags.push_back(std::sqrt(c[0]*c[0]+c[1]*c[1]+c[2]*c[2]));
        }
        if (!mags.empty()) {
            std::nth_element(mags.begin(), mags.begin()+mags.size()/2, mags.end());
            tau = std::max(0.05 * mags[mags.size()/2], 1e-12);
        }
    }

    // Per-collapse volume-change cap for Phase 2 (flap erase). A fold belongs
    // to a junk membrane only when it bounds *provably* zero volume — the two
    // sheets are coincident. A real folded-but-solid feature, however thin,
    // moves a measurable volume when collapsed. The cap is therefore set
    // near-exact-zero (1e-9·total): Phase 2 erases only provable junk and can
    // never damage real shape. The deliberate trade — a flap whose gap is not
    // quite zero is left intact (a leftover internal flap is still watertight
    // & manifold; a dented surface is a visible regression).
    double total_v6 = 0.0;
    for (const auto& f : F) total_v6 += tet6(V[f[0]], V[f[1]], V[f[2]]);
    const double vol_eps6 = std::max(std::abs(total_v6) * 1e-9,
                                     bbox_diag*bbox_diag*bbox_diag * 1e-15);

    struct Eval {
        std::vector<uint32_t> dead;
        std::vector<Tri>      born;
        int    nm_delta = 0, bnd_delta = 0, removed = 0;
        double vol6_delta = 0.0;
    };
    auto evaluate = [&](uint32_t u, uint32_t w, Eval& ev) -> bool {
        std::unordered_map<uint64_t, int> delta;
        auto it = v2f.find(w);
        if (it == v2f.end()) return false;
        for (uint32_t fi : it->second) {
            if (!alive[fi]) continue;
            const Tri& f = F[fi];
            if (f[0]!=w && f[1]!=w && f[2]!=w) continue;   // stale entry
            ev.dead.push_back(fi);
            ev.vol6_delta -= tet6(V[f[0]], V[f[1]], V[f[2]]);
            for (int k = 0; k < 3; ++k) delta[ekey(f[k], f[(k+1)%3])]--;
            Tri g{ f[0]==w?u:f[0], f[1]==w?u:f[1], f[2]==w?u:f[2] };
            if (g[0]==g[1] || g[1]==g[2] || g[0]==g[2]) { ++ev.removed; continue; }
            Vec3 oc = tri_cross(V[f[0]], V[f[1]], V[f[2]]);
            Vec3 nc = tri_cross(V[g[0]], V[g[1]], V[g[2]]);
            double oL = std::sqrt(oc[0]*oc[0]+oc[1]*oc[1]+oc[2]*oc[2]);
            double nL = std::sqrt(nc[0]*nc[0]+nc[1]*nc[1]+nc[2]*nc[2]);
            double d  = oc[0]*nc[0]+oc[1]*nc[1]+oc[2]*nc[2];
            if (d < 0.0 && oL > tau && nL > tau) return false;  // real flip
            ev.born.push_back(g);
            ev.vol6_delta += tet6(V[g[0]], V[g[1]], V[g[2]]);
            for (int k = 0; k < 3; ++k) delta[ekey(g[k], g[(k+1)%3])]++;
        }
        if (ev.dead.empty()) return false;
        for (const auto& kv : delta) {
            auto e = e2c.find(kv.first);
            int oldc = (e == e2c.end()) ? 0 : e->second;
            int newc = oldc + kv.second;
            ev.nm_delta  += (newc>2)  - (oldc>2);
            ev.bnd_delta += (newc==1) - (oldc==1);
        }
        return true;
    };
    auto commit = [&](const Eval& ev) {
        for (uint32_t fi : ev.dead) { sub_face(F[fi]); alive[fi] = 0; }
        for (const Tri& g : ev.born) {
            uint32_t gi = static_cast<uint32_t>(F.size());
            F.push_back(g);
            alive.push_back(1);
            add_face(g);
            for (int k = 0; k < 3; ++k) v2f[g[k]].push_back(gi);
        }
        nm  += ev.nm_delta;
        bnd += ev.bnd_delta;
    };

    int collapses = 0;

    // ---- Phase 1: erase non-manifold edges ----
    while (nm > 0 && collapses < max_collapses) {
        std::unordered_set<uint32_t> nmv;
        for (const auto& kv : e2c)
            if (kv.second > 2) {
                uint32_t a, b; edecode(kv.first, a, b);
                nmv.insert(a); nmv.insert(b);
            }
        std::vector<std::pair<double, uint64_t>> cand;
        for (const auto& kv : e2c) {
            uint32_t a, b; edecode(kv.first, a, b);
            if (!nmv.count(a) && !nmv.count(b)) continue;
            double L = vdist(V[a], V[b]);
            if (L <= max_collapse_len) cand.emplace_back(L, kv.first);
        }
        std::sort(cand.begin(), cand.end());

        bool did = false;
        for (int tier = 0; tier < 2 && !did; ++tier) {
            for (const auto& cd : cand) {
                uint32_t a, b; edecode(cd.second, a, b);
                for (int dir = 0; dir < 2; ++dir) {
                    Eval ev;
                    if (!evaluate(dir?b:a, dir?a:b, ev)) continue;
                    // Phase 1 relies on max_collapse_len + the whole-stage
                    // volume guard; a per-collapse cap here wrongly blocks
                    // legitimate NM collapses that nudge real wall geometry.
                    bool ok = (tier == 0)
                        ? (nm + ev.nm_delta < nm && bnd + ev.bnd_delta <= bnd)
                        : (ev.nm_delta <= 0 && ev.bnd_delta <= 0 && ev.removed > 0);
                    if (!ok) continue;
                    commit(ev);
                    ++collapses;
                    did = true;
                    break;
                }
                if (did) break;
            }
        }
        if (!did) break;
    }

    // ---- Phase 2: erase residual zero-volume flaps ----
    // collapse_nm Phase 1 removes the non-manifold edges but a doubled
    // membrane can survive as a now-manifold zero-volume flap — a fold where
    // two sheets lie back-to-back. Such a fold shows up as a manifold edge
    // whose two faces are near-anti-parallel. Collapsing those fold edges
    // erodes the flap; the per-collapse volume guard keeps real folded-but-
    // solid features (which enclose volume) untouched.
    if (nm == 0) {
        while (collapses < max_collapses) {
            // face normals (lazily, each pass — region is small)
            std::unordered_map<uint64_t, std::array<uint32_t,2>> ef;
            for (uint32_t fi = 0; fi < F.size(); ++fi) {
                if (!alive[fi]) continue;
                const Tri& f = F[fi];
                for (int k = 0; k < 3; ++k) {
                    uint64_t e = ekey(f[k], f[(k+1)%3]);
                    auto& slot = ef[e];
                    if (slot[0]==0 && slot[1]==0) slot = {fi+1u, 0u};
                    else if (slot[1]==0)          slot[1] = fi+1u;
                }
            }
            auto fnorm = [&](uint32_t fi) {
                const Tri& f = F[fi];
                Vec3 c = tri_cross(V[f[0]], V[f[1]], V[f[2]]);
                double L = std::sqrt(c[0]*c[0]+c[1]*c[1]+c[2]*c[2]);
                if (L < 1e-300) return Vec3{0,0,0};
                return Vec3{c[0]/L, c[1]/L, c[2]/L};
            };
            // fold edges: manifold edge whose 2 faces are near-anti-parallel.
            std::unordered_set<uint32_t> flapv;
            for (const auto& kv : ef) {
                if (kv.second[0]==0 || kv.second[1]==0) continue;
                uint32_t f1 = kv.second[0]-1, f2 = kv.second[1]-1;
                Vec3 n1 = fnorm(f1), n2 = fnorm(f2);
                double d = n1[0]*n2[0]+n1[1]*n2[1]+n1[2]*n2[2];
                if (d < -0.95) {
                    uint32_t a, b; edecode(kv.first, a, b);
                    flapv.insert(a); flapv.insert(b);
                }
            }
            if (flapv.empty()) break;

            std::vector<std::pair<double, uint64_t>> cand;
            for (const auto& kv : e2c) {
                uint32_t a, b; edecode(kv.first, a, b);
                if (!flapv.count(a) && !flapv.count(b)) continue;
                double L = vdist(V[a], V[b]);
                if (L <= max_collapse_len) cand.emplace_back(L, kv.first);
            }
            std::sort(cand.begin(), cand.end());

            bool did = false;
            for (const auto& cd : cand) {
                uint32_t a, b; edecode(cd.second, a, b);
                for (int dir = 0; dir < 2; ++dir) {
                    Eval ev;
                    if (!evaluate(dir?b:a, dir?a:b, ev)) continue;
                    if (std::abs(ev.vol6_delta) > vol_eps6) continue;
                    if (ev.nm_delta != 0 || ev.bnd_delta > 0 || ev.removed == 0)
                        continue;
                    commit(ev);
                    ++collapses;
                    did = true;
                    break;
                }
                if (did) break;
            }
            if (!did) break;
        }
    }

    Mesh out;
    out.vertices = mesh.vertices;
    out.faces.reserve(F.size());
    for (uint32_t fi = 0; fi < F.size(); ++fi)
        if (alive[fi]) out.faces.push_back(F[fi]);

    result.mesh      = std::move(out);
    result.nm_after  = static_cast<uint32_t>(nm);
    result.collapses = static_cast<uint32_t>(collapses);
    return result;
}

} // namespace meshseal::stages
