#include "thin.h"
#include "../internal/vec3.h"
#include <algorithm>
#include <cmath>
#include <unordered_map>
#include <vector>

namespace meshseal::stages {

namespace {

static uint64_t cell_key(int32_t cx, int32_t cy, int32_t cz) {
    // Pack three 21-bit signed ints into 64 bits
    auto pack = [](int32_t v) -> uint64_t {
        return static_cast<uint64_t>(static_cast<int64_t>(v) + (1 << 20)) & 0x1FFFFF;
    };
    return pack(cx) | (pack(cy) << 21) | (pack(cz) << 42);
}

} // anonymous namespace

ThinResult remove_thin_features(const Mesh& mesh,
                                double distance_threshold,
                                double cos_threshold) {
    using namespace internal;

    const uint32_t nf = static_cast<uint32_t>(mesh.faces.size());

    ThinResult result;
    result.faces_removed = 0;
    result.pairs_found   = 0;

    if (nf == 0) {
        result.mesh = mesh;
        return result;
    }

    // 1. Compute centroid and normal for each face
    std::vector<Vec3d> centroids(nf);
    std::vector<Vec3d> normals(nf);

    for (uint32_t i = 0; i < nf; ++i) {
        const auto& f = mesh.faces[i];
        const Vec3d& v0 = mesh.vertices[f[0]];
        const Vec3d& v1 = mesh.vertices[f[1]];
        const Vec3d& v2 = mesh.vertices[f[2]];

        centroids[i] = {
            (v0[0] + v1[0] + v2[0]) / 3.0,
            (v0[1] + v1[1] + v2[1]) / 3.0,
            (v0[2] + v1[2] + v2[2]) / 3.0
        };

        normals[i] = normalize(cross(sub(v1, v0), sub(v2, v0)));
    }

    // 2. Build spatial hash grid.
    //
    // Cell size needs two properties:
    //   - large enough that the 27-cell neighbourhood always contains every
    //     other face within distance_threshold (so we never miss a pair)
    //   - large enough that each cell holds only O(few) faces on a normal
    //     surface mesh (otherwise the 27-cell scan is dense and slow)
    //
    // The fixed 1e-3 default cell size used to violate the second property
    // on huge meshes: a 2 M-face mesh with mm-scale bbox has many faces
    // per mm³, so each 1 mm cell holds dozens of faces and the 27-cell
    // scan becomes O(face²)-ish locally. Adapt cell size to a multiple of
    // the mean edge length when that's larger than the requested
    // threshold: keeps the candidate set per cell small while still
    // catching every pair within distance_threshold (since 2 * mean_edge
    // already exceeds it on a uniformly-meshed surface).
    double mean_edge2 = 0.0;
    {
        const uint32_t sample_n = std::min<uint32_t>(nf, 4096);
        for (uint32_t i = 0; i < sample_n; ++i) {
            const auto& f = mesh.faces[i];
            const Vec3d& a = mesh.vertices[f[0]];
            const Vec3d& b = mesh.vertices[f[1]];
            const Vec3d& c = mesh.vertices[f[2]];
            const double e1 = (b[0]-a[0])*(b[0]-a[0]) + (b[1]-a[1])*(b[1]-a[1]) + (b[2]-a[2])*(b[2]-a[2]);
            const double e2 = (c[0]-b[0])*(c[0]-b[0]) + (c[1]-b[1])*(c[1]-b[1]) + (c[2]-b[2])*(c[2]-b[2]);
            const double e3 = (a[0]-c[0])*(a[0]-c[0]) + (a[1]-c[1])*(a[1]-c[1]) + (a[2]-c[2])*(a[2]-c[2]);
            mean_edge2 += (e1 + e2 + e3) / 3.0;
        }
        mean_edge2 /= sample_n;
    }
    const double mean_edge = std::sqrt(mean_edge2);
    const double thr = distance_threshold > 0.0 ? distance_threshold : 1e-3;
    // We need cell_size >= thr (so any pair within thr is in 27-neighbour
    // set) AND cell_size >= ~mean_edge (so cell occupancy stays small).
    const double cell_size = std::max(thr, mean_edge);

    std::unordered_map<uint64_t, std::vector<uint32_t>> grid;
    grid.reserve(nf);

    for (uint32_t i = 0; i < nf; ++i) {
        const auto& c = centroids[i];
        int32_t cx = static_cast<int32_t>(std::floor(c[0] / cell_size));
        int32_t cy = static_cast<int32_t>(std::floor(c[1] / cell_size));
        int32_t cz = static_cast<int32_t>(std::floor(c[2] / cell_size));
        grid[cell_key(cx, cy, cz)].push_back(i);
    }

    // 3. Find thin pairs
    std::vector<bool> marked(nf, false);
    const double dist2 = distance_threshold * distance_threshold;

    for (uint32_t i = 0; i < nf; ++i) {
        const auto& ci = centroids[i];
        int32_t bx = static_cast<int32_t>(std::floor(ci[0] / cell_size));
        int32_t by = static_cast<int32_t>(std::floor(ci[1] / cell_size));
        int32_t bz = static_cast<int32_t>(std::floor(ci[2] / cell_size));

        // Check all 27 neighboring cells (including own cell)
        for (int32_t dx = -1; dx <= 1; ++dx) {
            for (int32_t dy = -1; dy <= 1; ++dy) {
                for (int32_t dz = -1; dz <= 1; ++dz) {
                    auto it = grid.find(cell_key(bx + dx, by + dy, bz + dz));
                    if (it == grid.end()) continue;

                    for (uint32_t j : it->second) {
                        if (j <= i) continue; // only process each pair once (i < j)

                        const auto& cj = centroids[j];
                        double ddx = ci[0] - cj[0];
                        double ddy = ci[1] - cj[1];
                        double ddz = ci[2] - cj[2];
                        double d2 = ddx * ddx + ddy * ddy + ddz * ddz;

                        if (d2 > dist2) continue;

                        double dp = dot(normals[i], normals[j]);
                        if (dp < cos_threshold) {
                            marked[i] = true;
                            marked[j] = true;
                            ++result.pairs_found;
                        }
                    }
                }
            }
        }
    }

    // 4. Build output mesh excluding marked faces, compact vertices
    // Count how many faces survive
    std::vector<uint32_t> kept_faces;
    kept_faces.reserve(nf);
    for (uint32_t i = 0; i < nf; ++i) {
        if (!marked[i]) kept_faces.push_back(i);
    }
    result.faces_removed = nf - static_cast<uint32_t>(kept_faces.size());

    // Remap vertices: only include vertices referenced by surviving faces
    const uint32_t nv = static_cast<uint32_t>(mesh.vertices.size());
    std::vector<uint32_t> vert_remap(nv, UINT32_MAX);
    uint32_t new_idx = 0;

    for (uint32_t fi : kept_faces) {
        const auto& f = mesh.faces[fi];
        for (int k = 0; k < 3; ++k) {
            if (vert_remap[f[k]] == UINT32_MAX) {
                vert_remap[f[k]] = new_idx++;
            }
        }
    }

    Mesh out;
    out.vertices.resize(new_idx);
    for (uint32_t v = 0; v < nv; ++v) {
        if (vert_remap[v] != UINT32_MAX) {
            out.vertices[vert_remap[v]] = mesh.vertices[v];
        }
    }

    out.faces.reserve(kept_faces.size());
    for (uint32_t fi : kept_faces) {
        const auto& f = mesh.faces[fi];
        out.faces.push_back({vert_remap[f[0]], vert_remap[f[1]], vert_remap[f[2]]});
    }

    result.mesh = std::move(out);
    return result;
}

} // namespace meshseal::stages

