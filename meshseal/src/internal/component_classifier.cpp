#define _USE_MATH_DEFINES
#include <cmath>
#include "component_classifier.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <queue>
#include <unordered_map>
#include <unordered_set>

namespace meshseal::internal {

namespace {

// -------------------------------------------------------------------------
// Closed-form eigenvalues for a 3×3 symmetric matrix using the trigonometric
// (Cardano) method.  Returns eigenvalues in ascending order: e[0] ≤ e[1] ≤ e[2].
// Reference: Smith (1961) / Wikipedia "Eigenvalue algorithm for symmetric 3×3 matrices".
// -------------------------------------------------------------------------
static std::array<double, 3> sym3_eigenvalues(
    double a00, double a01, double a02,
    double a11, double a12,
    double a22)
{
    // Off-diagonal squared sum
    const double p1 = a01 * a01 + a02 * a02 + a12 * a12;
    if (p1 < 1e-30) {
        // A is diagonal (or numerically so)
        std::array<double, 3> e = {a00, a11, a22};
        std::sort(e.begin(), e.end());
        return e;
    }

    const double q   = (a00 + a11 + a22) / 3.0;
    const double p2  = (a00 - q) * (a00 - q)
                     + (a11 - q) * (a11 - q)
                     + (a22 - q) * (a22 - q)
                     + 2.0 * p1;
    const double p   = std::sqrt(p2 / 6.0);

    // Normalised matrix B = (A - q*I) / p
    const double b00 = (a00 - q) / p;
    const double b01 = a01 / p;
    const double b02 = a02 / p;
    const double b11 = (a11 - q) / p;
    const double b12 = a12 / p;
    const double b22 = (a22 - q) / p;

    // det(B) / 2
    const double r = 0.5 * (b00 * (b11 * b22 - b12 * b12)
                           - b01 * (b01 * b22 - b12 * b02)
                           + b02 * (b01 * b12 - b11 * b02));

    // Clamp r to [-1, 1] to guard against floating-point rounding errors in acos
    double phi;
    if      (r <= -1.0) phi = M_PI / 3.0;
    else if (r >=  1.0) phi = 0.0;
    else                phi = std::acos(r) / 3.0;

    std::array<double, 3> e;
    e[2] = q + 2.0 * p * std::cos(phi);                    // largest
    e[0] = q + 2.0 * p * std::cos(phi + 2.0 * M_PI / 3.0); // smallest
    e[1] = 3.0 * q - e[0] - e[2];                          // middle

    std::sort(e.begin(), e.end());
    return e;
}

// Encode an undirected edge as a single uint64_t key
static inline uint64_t edge_key(uint32_t a, uint32_t b) {
    if (a > b) std::swap(a, b);
    return (static_cast<uint64_t>(a) << 32) | static_cast<uint64_t>(b);
}

} // anonymous namespace

// -------------------------------------------------------------------------
// classify_components
// -------------------------------------------------------------------------
std::vector<ComponentInfo> classify_components(
    const Mesh& mesh,
    double planarity_thresh,
    double open_ratio_thresh)
{
    const uint32_t nf = static_cast<uint32_t>(mesh.faces.size());
    if (nf == 0) return {};

    // ------------------------------------------------------------------
    // Step 1: Build edge → face list and face adjacency.
    // ------------------------------------------------------------------
    // edge_valence: count faces per undirected edge (capped at 3 — we only need
    // to distinguish 1 (boundary) from ≥2 (interior/non-manifold)).
    std::unordered_map<uint64_t, uint8_t>               edge_valence;
    std::unordered_map<uint64_t, std::vector<uint32_t>> edge_to_faces;
    edge_valence.reserve(nf * 3);
    edge_to_faces.reserve(nf * 3);

    for (uint32_t fi = 0; fi < nf; ++fi) {
        const auto& f = mesh.faces[fi];
        for (int k = 0; k < 3; ++k) {
            const uint64_t ek = edge_key(f[k], f[(k + 1) % 3]);
            auto& val = edge_valence[ek];
            if (val < 3) ++val;
            edge_to_faces[ek].push_back(fi);
        }
    }

    // Face adjacency list (any shared edge connects two faces).
    std::vector<std::vector<uint32_t>> adj(nf);
    for (auto& [ek, inc] : edge_to_faces) {
        for (size_t i = 0; i < inc.size(); ++i)
            for (size_t j = i + 1; j < inc.size(); ++j) {
                adj[inc[i]].push_back(inc[j]);
                adj[inc[j]].push_back(inc[i]);
            }
    }

    // ------------------------------------------------------------------
    // Step 2: BFS to label connected components.
    // ------------------------------------------------------------------
    std::vector<int32_t> comp_label(nf, -1);
    int32_t num_components = 0;
    {
        std::queue<uint32_t> q;
        for (uint32_t seed = 0; seed < nf; ++seed) {
            if (comp_label[seed] != -1) continue;
            const int32_t c = num_components++;
            comp_label[seed] = c;
            q.push(seed);
            while (!q.empty()) {
                uint32_t cur = q.front(); q.pop();
                for (uint32_t nb : adj[cur])
                    if (comp_label[nb] == -1) { comp_label[nb] = c; q.push(nb); }
            }
        }
    }

    // ------------------------------------------------------------------
    // Step 3: Per-component statistics.
    //   - face_count
    //   - per-edge boundary count  → open_ratio = boundary_edges / total_edges
    //   - vertex positions        → PCA covariance → planarity_ratio
    // ------------------------------------------------------------------
    std::vector<ComponentInfo> infos(static_cast<size_t>(num_components));
    for (auto& info : infos) {
        info.face_count     = 0;
        info.open_ratio     = 0.0;
        info.planarity_ratio = 0.0;
        info.cls            = ComponentClass::NO_BOUNDARY; // default; set below
        info.reconstruction_attempted = false;
    }

    // Accumulate face indices per component
    for (uint32_t fi = 0; fi < nf; ++fi)
        infos[static_cast<size_t>(comp_label[fi])].face_indices.push_back(fi);

    // Fill face_count
    for (auto& info : infos)
        info.face_count = static_cast<uint32_t>(info.face_indices.size());

    // Per-component edge statistics
    // We count total_edge_slots and boundary_edge_slots per component.
    // Each face contributes 3 edge slots; a boundary edge (valence 1) is one slot
    // to one component; an interior edge (valence 2) contributes one slot to each
    // of its two incident faces' components (same component in a manifold mesh).
    // We track per-component (total_edges, boundary_edges) in unique-edge counts.
    {
        // For each component: track edges we have already visited
        std::vector<std::unordered_set<uint64_t>> comp_total_edges(
            static_cast<size_t>(num_components));
        std::vector<uint32_t> comp_boundary(static_cast<size_t>(num_components), 0u);

        for (uint32_t fi = 0; fi < nf; ++fi) {
            const int32_t c = comp_label[fi];
            const auto& f   = mesh.faces[fi];
            for (int k = 0; k < 3; ++k) {
                const uint64_t ek  = edge_key(f[k], f[(k + 1) % 3]);
                // Only count each undirected edge once per component
                if (comp_total_edges[static_cast<size_t>(c)].insert(ek).second) {
                    // First time this component sees this edge
                    if (edge_valence.at(ek) == 1u)
                        ++comp_boundary[static_cast<size_t>(c)];
                }
            }
        }

        for (int32_t c = 0; c < num_components; ++c) {
            const size_t total = comp_total_edges[static_cast<size_t>(c)].size();
            const uint32_t bnd = comp_boundary[static_cast<size_t>(c)];
            infos[static_cast<size_t>(c)].open_ratio =
                (total == 0) ? 0.0 : static_cast<double>(bnd) / static_cast<double>(total);
        }
    }

    // Per-component PCA planarity
    for (int32_t c = 0; c < num_components; ++c) {
        auto& info = infos[static_cast<size_t>(c)];

        // Collect unique vertices for this component
        std::unordered_set<uint32_t> vert_set;
        for (uint32_t fi : info.face_indices)
            for (int k = 0; k < 3; ++k)
                vert_set.insert(mesh.faces[fi][k]);

        const size_t nv = vert_set.size();
        if (nv < 2) { info.planarity_ratio = 0.0; continue; }

        // Compute mean
        double mx = 0, my = 0, mz = 0;
        for (uint32_t vi : vert_set) {
            mx += mesh.vertices[vi][0];
            my += mesh.vertices[vi][1];
            mz += mesh.vertices[vi][2];
        }
        const double inv_n = 1.0 / static_cast<double>(nv);
        mx *= inv_n; my *= inv_n; mz *= inv_n;

        // Accumulate 3×3 covariance (symmetric, 6 unique entries)
        double c00=0, c01=0, c02=0, c11=0, c12=0, c22=0;
        for (uint32_t vi : vert_set) {
            const double dx = mesh.vertices[vi][0] - mx;
            const double dy = mesh.vertices[vi][1] - my;
            const double dz = mesh.vertices[vi][2] - mz;
            c00 += dx * dx; c01 += dx * dy; c02 += dx * dz;
            c11 += dy * dy; c12 += dy * dz;
            c22 += dz * dz;
        }
        c00 *= inv_n; c01 *= inv_n; c02 *= inv_n;
        c11 *= inv_n; c12 *= inv_n; c22 *= inv_n;

        auto evals = sym3_eigenvalues(c00, c01, c02, c11, c12, c22);
        // evals[0] ≤ evals[1] ≤ evals[2]
        // Guard against degenerate (all-same-point) mesh
        info.planarity_ratio = (evals[2] < 1e-30)
            ? 0.0
            : evals[0] / evals[2];
    }

    // ------------------------------------------------------------------
    // Step 4: Apply routing rules (first match wins).
    // ------------------------------------------------------------------
    for (auto& info : infos) {
        if (info.open_ratio >= 1.0 - 1e-9) {
            // all edges are boundary edges → geometrically disconnected fragments
            info.cls = ComponentClass::SOUP;
        } else if (info.planarity_ratio < planarity_thresh &&
                   info.open_ratio > open_ratio_thresh) {
            info.cls = ComponentClass::SOUP;
        } else if (info.open_ratio > 1e-9) {
            info.cls = ComponentClass::OPEN;
        } else {
            info.cls = ComponentClass::NO_BOUNDARY;
        }
    }

    return infos;
}

// -------------------------------------------------------------------------
// extract_component
// -------------------------------------------------------------------------
Mesh extract_component(const Mesh& src, const std::vector<uint32_t>& face_indices)
{
    // Build remapping: old vertex index → new compact index
    std::unordered_map<uint32_t, uint32_t> remap;
    remap.reserve(face_indices.size() * 3);

    Mesh out;
    out.faces.reserve(face_indices.size());

    for (uint32_t fi : face_indices) {
        const auto& f = src.faces[fi];
        std::array<uint32_t, 3> nf{};
        for (int k = 0; k < 3; ++k) {
            auto it = remap.find(f[k]);
            if (it == remap.end()) {
                const uint32_t new_idx = static_cast<uint32_t>(out.vertices.size());
                remap[f[k]] = new_idx;
                out.vertices.push_back(src.vertices[f[k]]);
                nf[k] = new_idx;
            } else {
                nf[k] = it->second;
            }
        }
        out.faces.push_back(nf);
    }

    return out;
}

// -------------------------------------------------------------------------
// merge_meshes
// -------------------------------------------------------------------------
Mesh merge_meshes(const std::vector<Mesh>& meshes)
{
    Mesh out;
    size_t total_verts = 0, total_faces = 0;
    for (const auto& m : meshes) {
        total_verts += m.vertices.size();
        total_faces += m.faces.size();
    }
    out.vertices.reserve(total_verts);
    out.faces.reserve(total_faces);

    for (const auto& m : meshes) {
        const uint32_t v_offset = static_cast<uint32_t>(out.vertices.size());
        for (const auto& v : m.vertices)
            out.vertices.push_back(v);
        for (const auto& f : m.faces)
            out.faces.push_back({f[0] + v_offset, f[1] + v_offset, f[2] + v_offset});
    }

    return out;
}

} // namespace meshseal::internal
