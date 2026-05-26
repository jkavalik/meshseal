#include "orient_wn.h"
#include "../internal/winding_number.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <queue>
#include <unordered_map>
#include <vector>

namespace meshseal::stages {

namespace {

inline uint64_t ekey(uint32_t a, uint32_t b) {
    if (a > b) std::swap(a, b);
    return (static_cast<uint64_t>(a) << 32) | static_cast<uint64_t>(b);
}

// Per-face normal (not normalised — caller normalises if needed).
inline std::array<double, 3> face_normal_raw(const Mesh& m, uint32_t fi) {
    const auto& f = m.faces[fi];
    const auto& a = m.vertices[f[0]];
    const auto& b = m.vertices[f[1]];
    const auto& c = m.vertices[f[2]];
    return {
        (b[1]-a[1]) * (c[2]-a[2]) - (b[2]-a[2]) * (c[1]-a[1]),
        (b[2]-a[2]) * (c[0]-a[0]) - (b[0]-a[0]) * (c[2]-a[2]),
        (b[0]-a[0]) * (c[1]-a[1]) - (b[1]-a[1]) * (c[0]-a[0]),
    };
}

inline double face_area_doubled(const std::array<double, 3>& n) {
    return std::sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);
}

inline std::array<double, 3> face_centroid(const Mesh& m, uint32_t fi) {
    const auto& f = m.faces[fi];
    const auto& a = m.vertices[f[0]];
    const auto& b = m.vertices[f[1]];
    const auto& c = m.vertices[f[2]];
    return {(a[0]+b[0]+c[0]) / 3.0,
            (a[1]+b[1]+c[1]) / 3.0,
            (a[2]+b[2]+c[2]) / 3.0};
}

// Count manifold edges where the two incident faces have near-anti-
// parallel normals (dot < -0.95). Same definition as the diagnostic.
uint32_t count_antipar_manifold_pairs(const Mesh& m) {
    const uint32_t nf = static_cast<uint32_t>(m.faces.size());
    std::unordered_map<uint64_t, std::vector<uint32_t>> edge_to_faces;
    edge_to_faces.reserve(nf * 3);
    for (uint32_t fi = 0; fi < nf; ++fi) {
        const auto& f = m.faces[fi];
        for (int k = 0; k < 3; ++k) {
            const uint32_t u = f[k], v = f[(k+1) % 3];
            edge_to_faces[ekey(u, v)].push_back(fi);
        }
    }
    // Pre-compute normalised normals.
    std::vector<std::array<double, 3>> n(nf);
    for (uint32_t fi = 0; fi < nf; ++fi) {
        auto raw = face_normal_raw(m, fi);
        double L = face_area_doubled(raw);
        if (L < 1e-30) {
            n[fi] = {0.0, 0.0, 0.0};
        } else {
            n[fi] = {raw[0]/L, raw[1]/L, raw[2]/L};
        }
    }
    uint32_t antipar = 0;
    for (const auto& kv : edge_to_faces) {
        if (kv.second.size() != 2) continue;
        const uint32_t a = kv.second[0], b = kv.second[1];
        const double dot = n[a][0]*n[b][0] + n[a][1]*n[b][1] + n[a][2]*n[b][2];
        if (dot < -0.95) ++antipar;
    }
    return antipar;
}

} // anonymous namespace

OrientWnResult orient_by_winding_number(const Mesh& mesh,
                                        uint32_t min_antipar_pairs,
                                        size_t   max_faces)
{
    OrientWnResult result;
    result.mesh = mesh;
    const uint32_t nf = static_cast<uint32_t>(mesh.faces.size());
    if (nf == 0 || nf > max_faces) return result;

    // Count NM edges (will be re-counted post-flip as a guard).
    auto count_nm_edges = [](const Mesh& m) {
        std::unordered_map<uint64_t, int> ec;
        for (const auto& f : m.faces) {
            for (int k = 0; k < 3; ++k) {
                ec[ekey(f[k], f[(k+1) % 3])] += 1;
            }
        }
        uint32_t nm = 0;
        for (const auto& kv : ec) if (kv.second > 2) ++nm;
        return nm;
    };
    const uint32_t nm_pre = count_nm_edges(mesh);
    // Skip when the mesh still has NM edges. orient_wn is designed
    // for the post-cleanup state. On meshes with NM, the patch BFS
    // can't cleanly separate orientation regions and the relative-GWN
    // test gives ambiguous results near the NM boundaries — the
    // resulting flips can subtly degrade downstream stages even when
    // local antipar checks pass. Concrete observation (kuzely.stl,
    // 2026-05-25): orient_wn in a recursive carve_refill call flipped
    // 3 faces, antipar 1337->1336 strictly decreased, but downstream
    // pipeline ended with nm=1 (was CLEAN before).
    if (nm_pre > 0) {
        result.antipar_before = count_antipar_manifold_pairs(mesh);
        result.antipar_after  = result.antipar_before;
        return result;
    }

    // Gate: count antipar pairs. Skip if below threshold.
    const uint32_t antipar_pre = count_antipar_manifold_pairs(mesh);
    result.antipar_before = antipar_pre;
    if (antipar_pre < min_antipar_pairs) {
        result.antipar_after = antipar_pre;
        return result;
    }

    // Build undirected edge → faces map; mark non-manifold edges.
    std::unordered_map<uint64_t, std::vector<uint32_t>> edge_to_faces;
    edge_to_faces.reserve(nf * 3);
    for (uint32_t fi = 0; fi < nf; ++fi) {
        const auto& f = mesh.faces[fi];
        for (int k = 0; k < 3; ++k) {
            const uint32_t u = f[k], v = f[(k+1) % 3];
            edge_to_faces[ekey(u, v)].push_back(fi);
        }
    }

    // Pre-compute per-face raw normal and area.
    std::vector<std::array<double, 3>> nraw(nf);
    std::vector<double> area(nf);
    for (uint32_t fi = 0; fi < nf; ++fi) {
        nraw[fi] = face_normal_raw(mesh, fi);
        area[fi] = 0.5 * face_area_doubled(nraw[fi]);
    }

    // Build face_neighbors for ORIENTATION-PATCH BFS: connect faces
    // sharing an edge with OPPOSITE winding (= parallel normals).
    // Faces sharing an edge with SAME winding (anti-parallel normals)
    // are NOT connected — they belong to different patches.
    struct FaceEdge { uint32_t neighbor; uint32_t u, v; };
    std::vector<std::vector<FaceEdge>> face_neighbors(nf);
    for (uint32_t fi = 0; fi < nf; ++fi) {
        const auto& f = mesh.faces[fi];
        for (int k = 0; k < 3; ++k) {
            const uint32_t u = f[k], v = f[(k+1) % 3];
            const uint64_t key = ekey(u, v);
            const auto& shared = edge_to_faces.at(key);
            if (shared.size() != 2) continue;  // skip NM edges
            for (uint32_t nb : shared) {
                if (nb != fi) face_neighbors[fi].push_back({nb, u, v});
            }
        }
    }

    // BFS to find orientation patches.
    std::vector<int> patch(nf, -1);
    std::vector<std::vector<uint32_t>> patch_faces;
    for (uint32_t seed = 0; seed < nf; ++seed) {
        if (patch[seed] != -1) continue;
        const int pid = static_cast<int>(patch_faces.size());
        patch_faces.emplace_back();
        patch[seed] = pid;
        std::queue<uint32_t> q;
        q.push(seed);
        while (!q.empty()) {
            const uint32_t cur = q.front(); q.pop();
            patch_faces[pid].push_back(cur);
            for (const auto& fe : face_neighbors[cur]) {
                const uint32_t nb = fe.neighbor;
                if (patch[nb] != -1) continue;
                // Connected if winding is OPPOSITE (parallel normals).
                const auto& fc = mesh.faces[cur];
                bool cur_has_uv = false;
                for (int e = 0; e < 3; ++e) {
                    if (fc[e] == fe.u && fc[(e + 1) % 3] == fe.v) {
                        cur_has_uv = true; break;
                    }
                }
                const uint32_t real_u = cur_has_uv ? fe.u : fe.v;
                const uint32_t real_v = cur_has_uv ? fe.v : fe.u;
                const auto& fn = mesh.faces[nb];
                bool nb_has_vu = false;
                for (int e = 0; e < 3; ++e) {
                    if (fn[e] == real_v && fn[(e + 1) % 3] == real_u) {
                        nb_has_vu = true; break;
                    }
                }
                if (nb_has_vu) {
                    patch[nb] = pid;
                    q.push(nb);
                }
            }
        }
    }

    // bbox for ε scaling.
    std::array<double, 3> lo{mesh.vertices[0][0], mesh.vertices[0][1], mesh.vertices[0][2]};
    std::array<double, 3> hi = lo;
    for (const auto& v : mesh.vertices)
        for (int k = 0; k < 3; ++k) {
            if (v[k] < lo[k]) lo[k] = v[k];
            if (v[k] > hi[k]) hi[k] = v[k];
        }
    const double dx = hi[0]-lo[0], dy = hi[1]-lo[1], dz = hi[2]-lo[2];
    const double bbox_diag = std::sqrt(dx*dx + dy*dy + dz*dz);
    const double eps = bbox_diag * 1e-5;

    // Per-face flip-decision: use the RELATIVE GWN comparison. For each
    // face F sample two points: p_out = centroid + eps·N (where the
    // normal SAYS outside is) and p_in = centroid - eps·N (where the
    // normal SAYS inside is). For a correctly-oriented face,
    // GWN(p_out) < GWN(p_in). For an inverted face, GWN(p_out) >
    // GWN(p_in). The DIFFERENCE is meaningful even when both absolute
    // GWNs are fractional (membrane interior case).
    //
    // For large patches, sub-sample faces (up to kMaxSamplesPerPatch
    // random ones) and majority-vote. For small patches, evaluate all
    // faces. Adopt the flip iff the majority disagrees with the patch's
    // current orientation.
    constexpr size_t kMaxSamplesPerPatch = 30;
    std::vector<bool> flip_patch(patch_faces.size(), false);
    auto face_signed_vol = [&](const Mesh& m, uint32_t fi) {
        const auto& f = m.faces[fi];
        const auto& v0 = m.vertices[f[0]];
        const auto& v1 = m.vertices[f[1]];
        const auto& v2 = m.vertices[f[2]];
        return (v0[0] * (v1[1] * v2[2] - v1[2] * v2[1])
              - v0[1] * (v1[0] * v2[2] - v1[2] * v2[0])
              + v0[2] * (v1[0] * v2[1] - v1[1] * v2[0])) / 6.0;
    };
    // Pre-flip total signed volume — Guard 2 below rejects any flip that
    // inverts this sign (wholesale-inverted result, the Bee_v3 vol[3]
    // failure mode that motivated the relative-GWN refinement).
    double vol_pre = 0.0;
    for (uint32_t fi = 0; fi < nf; ++fi) vol_pre += face_signed_vol(mesh, fi);

    // Simple LCG for deterministic sub-sampling.
    auto sample_step = [](uint64_t& s) {
        s = s * 6364136223846793005ull + 1442695040888963407ull;
        return s;
    };
    for (size_t pid = 0; pid < patch_faces.size(); ++pid) {
        const auto& fs = patch_faces[pid];
        if (fs.empty()) continue;
        // Decide how many faces to sample.
        std::vector<uint32_t> samples;
        if (fs.size() <= kMaxSamplesPerPatch) {
            samples = fs;
        } else {
            samples.reserve(kMaxSamplesPerPatch);
            uint64_t s = static_cast<uint64_t>(pid) * 0x9E3779B97F4A7C15ull + 1;
            for (size_t k = 0; k < kMaxSamplesPerPatch; ++k) {
                samples.push_back(fs[sample_step(s) % fs.size()]);
            }
        }
        int votes_flip = 0;
        int votes_keep = 0;
        for (uint32_t fi : samples) {
            if (area[fi] < 1e-30) continue;  // skip slivers
            auto c = face_centroid(mesh, fi);
            auto n = nraw[fi];
            const double L = face_area_doubled(n);
            n[0] /= L; n[1] /= L; n[2] /= L;
            std::array<double, 3> p_out{
                c[0] + eps * n[0], c[1] + eps * n[1], c[2] + eps * n[2]};
            std::array<double, 3> p_in {
                c[0] - eps * n[0], c[1] - eps * n[1], c[2] - eps * n[2]};
            const double wn_out = internal::generalized_winding_number(mesh, p_out);
            const double wn_in  = internal::generalized_winding_number(mesh, p_in);
            // If wn_out > wn_in: the "outside" point is more inside than
            // the "inside" point — face is inverted, vote flip. Tolerance
            // 1e-3 to ignore near-tie cases (likely thin features or
            // numerical noise).
            if (wn_out > wn_in + 1e-3) ++votes_flip;
            else if (wn_in > wn_out + 1e-3) ++votes_keep;
        }
        // Need a clear majority to flip: at least 60% of the deciding
        // votes, AND at least 3 votes for flip (avoid coincidence on
        // tiny patches).
        const int total = votes_flip + votes_keep;
        if (total >= 3 && votes_flip * 5 >= total * 3) {
            flip_patch[pid] = true;
        }
    }

    // Apply flips.
    auto out_mesh = mesh;
    uint32_t flipped = 0;
    for (size_t pid = 0; pid < patch_faces.size(); ++pid) {
        if (!flip_patch[pid]) continue;
        for (uint32_t fi : patch_faces[pid]) {
            std::swap(out_mesh.faces[fi][1], out_mesh.faces[fi][2]);
            ++flipped;
        }
    }
    if (flipped == 0) {
        result.antipar_after = antipar_pre;
        return result;
    }

    // Guard 1: antipar count must strictly decrease.
    const uint32_t antipar_post = count_antipar_manifold_pairs(out_mesh);
    result.antipar_after = antipar_post;
    if (antipar_post >= antipar_pre) {
        // Revert.
        return result;
    }
    // Guard 1b: NM edge count must not increase. (Flipping a face shouldn't
    // change edge incidence counts in principle, but on inputs with
    // pre-existing NM edges or float-precision quirks, the topology can
    // become more fragile. Be strict.)
    const uint32_t nm_post = count_nm_edges(out_mesh);
    if (nm_post > nm_pre) {
        return result;
    }
    // Guard 2: total signed-volume sign must not flip (wholesale
    // inversion — the Bee_v3 vol[3] failure mode).
    double vol_post = 0.0;
    for (uint32_t fi = 0; fi < nf; ++fi) vol_post += face_signed_vol(out_mesh, fi);
    if (std::abs(vol_pre) > 1e-9 && vol_pre * vol_post < 0.0) {
        return result;  // revert: would invert the solid
    }
    result.mesh = std::move(out_mesh);
    result.applied = true;
    result.faces_flipped = flipped;
    return result;
}

} // namespace meshseal::stages
