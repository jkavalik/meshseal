#include "tjunction.h"
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

} // anonymous namespace

TJunctionResult split_tjunctions(const Mesh& mesh) {
    TJunctionResult result;
    result.mesh = mesh;

    const uint32_t nf = static_cast<uint32_t>(mesh.faces.size());
    if (nf == 0 || mesh.vertices.empty()) return result;

    const auto& V = mesh.vertices;

    // bbox diag → collinearity tolerance.
    Vec3 lo = V[0], hi = V[0];
    for (const auto& p : V)
        for (int k = 0; k < 3; ++k) {
            lo[k] = std::min(lo[k], p[k]);
            hi[k] = std::max(hi[k], p[k]);
        }
    double dx = hi[0]-lo[0], dy = hi[1]-lo[1], dz = hi[2]-lo[2];
    const double bbox_diag = std::sqrt(dx*dx + dy*dy + dz*dz);
    const double perp_tol = bbox_diag * 1e-6;

    // Edge → incident faces.
    std::unordered_map<uint64_t, std::vector<uint32_t>> e2f;
    e2f.reserve(nf * 3);
    for (uint32_t fi = 0; fi < nf; ++fi) {
        const auto& f = mesh.faces[fi];
        for (int k = 0; k < 3; ++k)
            e2f[ekey(f[k], f[(k+1)%3])].push_back(fi);
    }

    // Problem edges = boundary or non-manifold (count != 2); collect their
    // vertices as the candidate set for the collinearity scan.
    std::vector<uint64_t> problem_edges;
    std::unordered_set<uint32_t> cand_set;
    for (const auto& kv : e2f) {
        if (kv.second.size() == 2) continue;
        problem_edges.push_back(kv.first);
        uint32_t a = static_cast<uint32_t>(kv.first >> 32);
        uint32_t b = static_cast<uint32_t>(kv.first & 0xffffffffu);
        cand_set.insert(a);
        cand_set.insert(b);
    }
    if (problem_edges.empty()) return result;
    std::vector<uint32_t> cand(cand_set.begin(), cand_set.end());

    // For each problem edge, find vertices collinear on its interior.
    // splits[fi] : for face fi, the ordered insert-chain on one of its edges.
    struct FaceSplit {
        uint32_t a, b;                  // the edge being split (as found)
        std::vector<uint32_t> chain;    // a, v1, …, vk, b  (ordered along a→b)
    };
    std::unordered_map<uint32_t, FaceSplit> face_split;

    for (uint64_t ek : problem_edges) {
        uint32_t a = static_cast<uint32_t>(ek >> 32);
        uint32_t b = static_cast<uint32_t>(ek & 0xffffffffu);
        const Vec3& pa = V[a];
        const Vec3& pb = V[b];
        double d[3] = { pb[0]-pa[0], pb[1]-pa[1], pb[2]-pa[2] };
        double dd = d[0]*d[0] + d[1]*d[1] + d[2]*d[2];
        if (dd < 1e-300) continue;

        std::vector<std::pair<double,uint32_t>> on_edge;  // (t, vertex)
        for (uint32_t v : cand) {
            if (v == a || v == b) continue;
            const Vec3& pv = V[v];
            double e[3] = { pv[0]-pa[0], pv[1]-pa[1], pv[2]-pa[2] };
            double t = (e[0]*d[0] + e[1]*d[1] + e[2]*d[2]) / dd;
            if (t <= 1e-6 || t >= 1.0 - 1e-6) continue;   // not interior
            // perpendicular distance of v from the line
            double cx = d[1]*e[2] - d[2]*e[1];
            double cy = d[2]*e[0] - d[0]*e[2];
            double cz = d[0]*e[1] - d[1]*e[0];
            double perp = std::sqrt(cx*cx + cy*cy + cz*cz) / std::sqrt(dd);
            if (perp > perp_tol) continue;
            on_edge.emplace_back(t, v);
        }
        if (on_edge.empty()) continue;
        std::sort(on_edge.begin(), on_edge.end());

        std::vector<uint32_t> chain;
        chain.push_back(a);
        for (const auto& tv : on_edge) {
            if (chain.back() != tv.second) chain.push_back(tv.second);
        }
        chain.push_back(b);
        if (chain.size() < 3) continue;   // nothing inserted

        ++result.edges_split;
        for (uint32_t fi : e2f[ek]) {
            // Keep the first split found for a face (rare to have two).
            if (face_split.find(fi) == face_split.end())
                face_split[fi] = FaceSplit{ a, b, chain };
        }
    }

    if (face_split.empty()) { result.edges_split = 0; return result; }

    // Rebuild the face list, fanning split faces through their chain.
    Mesh out;
    out.vertices = mesh.vertices;
    out.faces.reserve(nf + face_split.size() * 3);
    for (uint32_t fi = 0; fi < nf; ++fi) {
        auto it = face_split.find(fi);
        if (it == face_split.end()) {
            out.faces.push_back(mesh.faces[fi]);
            continue;
        }
        const Tri& f = mesh.faces[fi];
        const FaceSplit& fs = it->second;
        // Locate the edge (fs.a, fs.b) within the face and its apex.
        int ia = -1;
        for (int k = 0; k < 3; ++k) if (f[k] == fs.a) ia = k;
        if (ia < 0) { out.faces.push_back(f); continue; }
        const bool forward = (f[(ia+1)%3] == fs.b);   // face runs a→b
        const bool backward = (f[(ia+2)%3] == fs.b);  // face runs b→a
        if (!forward && !backward) { out.faces.push_back(f); continue; }
        uint32_t apex = f[(ia+ (forward ? 2 : 1)) % 3];
        const auto& ch = fs.chain;   // a … b
        for (size_t i = 0; i + 1 < ch.size(); ++i) {
            if (forward) out.faces.push_back({ ch[i], ch[i+1], apex });
            else         out.faces.push_back({ ch[i+1], ch[i], apex });
            ++result.faces_split;
        }
    }

    result.mesh = std::move(out);
    return result;
}

} // namespace meshseal::stages
