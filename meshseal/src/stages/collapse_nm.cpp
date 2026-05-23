#include "collapse_nm.h"
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
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
    // nm_edges: the set of non-manifold edges, maintained incrementally so
    // each collapse iteration does NOT rescan the whole edge map. Without
    // this the Phase-1 loop is O(E) per collapse → O(E·collapses), which
    // is multi-minute on 400k-face meshes (profiled: kytka1).
    std::unordered_set<uint64_t> nm_edges;
    for (const auto& kv : e2c) {
        if (kv.second > 2) { ++nm; nm_edges.insert(kv.first); }
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
        std::unordered_map<uint64_t, int> delta;  // per-edge count change
    };
    auto evaluate = [&](uint32_t u, uint32_t w, Eval& ev) -> bool {
        auto& delta = ev.delta;
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
    // v2f compactor: commit appends to v2f but never removes dead/stale
    // entries. On long collapse runs v2f[hot-vertex] bloats with garbage,
    // making evaluate (which iterates v2f[w]) and gather_cand (which
    // iterates v2f[x] for x in nmv/flapv) slower with every collapse.
    // Compact when total v2f entries exceed 2× the live-entry high water.
    std::size_t v2f_entries = 0;
    for (const auto& kv : v2f) v2f_entries += kv.second.size();
    std::size_t live_v2f_target = v2f_entries;
    auto compact_v2f = [&]() {
        std::size_t new_total = 0;
        for (auto& kv : v2f) {
            const uint32_t x = kv.first;
            auto& vec = kv.second;
            auto wr = vec.begin();
            for (auto rd = vec.begin(); rd != vec.end(); ++rd) {
                const uint32_t fi = *rd;
                if (!alive[fi]) continue;
                const Tri& f = F[fi];
                if (f[0] != x && f[1] != x && f[2] != x) continue;
                *wr++ = fi;
            }
            vec.erase(wr, vec.end());
            new_total += vec.size();
        }
        v2f_entries = new_total;
        live_v2f_target = new_total;
    };

    auto commit = [&](const Eval& ev) {
        for (uint32_t fi : ev.dead) { sub_face(F[fi]); alive[fi] = 0; }
        for (const Tri& g : ev.born) {
            uint32_t gi = static_cast<uint32_t>(F.size());
            F.push_back(g);
            alive.push_back(1);
            add_face(g);
            for (int k = 0; k < 3; ++k) {
                v2f[g[k]].push_back(gi);
                ++v2f_entries;
            }
        }
        nm  += ev.nm_delta;
        bnd += ev.bnd_delta;
        // Keep nm_edges in sync — only the edges in ev.delta can have
        // changed count; everything else is untouched. O(|delta|).
        for (const auto& kv : ev.delta) {
            auto e = e2c.find(kv.first);
            const int c = (e == e2c.end()) ? 0 : e->second;
            if (c > 2) nm_edges.insert(kv.first);
            else       nm_edges.erase(kv.first);
        }
    };

    // Gather collapse candidates: short edges incident to any vertex in
    // `vset`, reached via the v2f adjacency (NOT a full e2c scan). Produces
    // exactly the edge set a full scan would — every e2c edge sits on an
    // alive face, which is in v2f of its vertices — and the (L, key) sort
    // makes the result order-independent of how edges are visited.
    auto gather_cand = [&](const std::unordered_set<uint32_t>& vset,
                           std::vector<std::pair<double, uint64_t>>& cand) {
        std::unordered_set<uint64_t> seen;
        for (uint32_t x : vset) {
            auto it = v2f.find(x);
            if (it == v2f.end()) continue;
            for (uint32_t fi : it->second) {
                if (!alive[fi]) continue;
                const Tri& f = F[fi];
                if (f[0]!=x && f[1]!=x && f[2]!=x) continue;   // stale
                for (int k = 0; k < 3; ++k) {
                    const uint64_t e = ekey(f[k], f[(k+1)%3]);
                    if (!seen.insert(e).second) continue;
                    uint32_t a, b; edecode(e, a, b);
                    if (!vset.count(a) && !vset.count(b)) continue;
                    const double L = vdist(V[a], V[b]);
                    if (L <= max_collapse_len) cand.emplace_back(L, e);
                }
            }
        }
        std::sort(cand.begin(), cand.end());
    };

    int collapses = 0;

    const bool _prof = std::getenv("MESHSEAL_PROFILE") != nullptr;
    auto _t0 = std::chrono::steady_clock::now();
    auto _elapsed = [&]() {
        return std::chrono::duration<double>(
            std::chrono::steady_clock::now() - _t0).count();
    };
    auto _phaselog = [&](const char* label, int colls, int nm_now) {
        if (_prof) std::fprintf(stderr, "[cnm] %-10s elapsed=%6.2fs collapses=%d nm=%d F=%zu\n",
                                label, _elapsed(), colls, nm_now, F.size());
    };
    _phaselog("ph1_start", 0, nm);

    // ---- Phase 1: erase non-manifold edges ----
    // Stall detector: the "neutral" tier (tier=1) accepts collapses that
    // hold nm but decimate the region — useful occasionally to escape a
    // local stall, but on hard real-world inputs (e.g. kytka1, 400k faces,
    // nm=189) it can spiral indefinitely with little or no nm progress.
    // Cap consecutive non-improving collapses at 50: that's plenty for the
    // neutral tier to break a genuine stall, beyond which the stall is
    // structural and we should accept the partial result.
    int ph1_iter = 0;
    int stall_streak = 0;
    constexpr int kStallCap = 50;
    while (nm > 0 && collapses < max_collapses) {
        if (_prof && (ph1_iter % 100 == 0))
            std::fprintf(stderr, "[cnm] ph1 iter=%d elapsed=%6.2fs collapses=%d nm=%d F=%zu nm_edges=%zu\n",
                         ph1_iter, _elapsed(), collapses, nm, F.size(), nm_edges.size());
        ++ph1_iter;
        const int nm_at_iter_start = nm;
        // NM vertices = endpoints of the incrementally-maintained nm_edges
        // set (no full e2c scan).
        std::unordered_set<uint32_t> nmv;
        for (uint64_t e : nm_edges) {
            uint32_t a, b; edecode(e, a, b);
            nmv.insert(a); nmv.insert(b);
        }
        std::vector<std::pair<double, uint64_t>> cand;
        gather_cand(nmv, cand);

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
        if (nm < nm_at_iter_start) stall_streak = 0;
        else if (++stall_streak >= kStallCap) {
            if (_prof) std::fprintf(stderr,
                "[cnm] ph1 stall break at iter=%d collapses=%d nm=%d (%d consecutive non-improving)\n",
                ph1_iter, collapses, nm, stall_streak);
            break;
        }
        // Compact v2f when bloat exceeds 2× the high-water-mark of live
        // entries; otherwise v2f[hot-vertex] grows unboundedly and slows
        // every subsequent evaluate / gather_cand call.
        if (v2f_entries > live_v2f_target * 2 && v2f_entries > 1024) {
            compact_v2f();
        }
    }

    _phaselog("ph1_done", collapses, nm);

    // ---- Phase 2: erase residual zero-volume flaps ----
    // collapse_nm Phase 1 removes the non-manifold edges but a doubled
    // membrane can survive as a now-manifold zero-volume flap — a fold where
    // two sheets lie back-to-back. Such a fold shows up as a manifold edge
    // whose two faces are near-anti-parallel. Collapsing those fold edges
    // erodes the flap; the per-collapse volume guard keeps real folded-but-
    // solid features (which enclose volume) untouched.
    if (nm == 0) {
        int ph2_iters = 0;
        while (collapses < max_collapses) {
            if (_prof && (ph2_iters % 100 == 0))
                std::fprintf(stderr, "[cnm] ph2 iter=%d elapsed=%6.2fs collapses=%d F=%zu\n",
                             ph2_iters, _elapsed(), collapses, F.size());
            ++ph2_iters;
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
            gather_cand(flapv, cand);

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
            if (v2f_entries > live_v2f_target * 2 && v2f_entries > 1024) {
                compact_v2f();
            }
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
