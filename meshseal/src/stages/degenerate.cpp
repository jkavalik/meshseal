#include "degenerate.h"
#include "../internal/vec3.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <unordered_set>
#include <vector>

namespace meshseal::stages {

Mesh remove_degenerate_faces(const Mesh& mesh, double area_threshold) {
    // Determine effective threshold:
    //   - If caller passed an explicit absolute value (anything other than
    //     the sentinel default 0.0), use it as-is.
    //   - Otherwise, derive a scale-aware threshold from the mesh bbox:
    //     1e-22 of bbox-diagonal² — well below any meaningful triangle
    //     area (incl. Liepa's narrow slivers, which are valid filling
    //     triangles even when their absolute area is small), but well
    //     above floating-point noise (~1e-32 × bbox² for double).
    //
    // Why this matters: an absolute threshold like 1e-14 worked for
    // unit-scale meshes but mis-filtered legitimate sliver triangles
    // produced by Liepa hole-fill on real-world meshes (mm-to-meter
    // scale). The scale-aware default keeps those slivers while still
    // catching truly collinear zero-area triangles.
    double effective = area_threshold;
    if (effective == 0.0) {
        double scale_sq = 1.0;
        if (!mesh.vertices.empty()) {
            double lo[3] = { std::numeric_limits<double>::max(),
                             std::numeric_limits<double>::max(),
                             std::numeric_limits<double>::max() };
            double hi[3] = { std::numeric_limits<double>::lowest(),
                             std::numeric_limits<double>::lowest(),
                             std::numeric_limits<double>::lowest() };
            for (const auto& v : mesh.vertices) {
                for (int k = 0; k < 3; ++k) {
                    lo[k] = std::min(lo[k], v[k]);
                    hi[k] = std::max(hi[k], v[k]);
                }
            }
            const double dx = hi[0] - lo[0];
            const double dy = hi[1] - lo[1];
            const double dz = hi[2] - lo[2];
            const double d2 = dx * dx + dy * dy + dz * dz;
            if (d2 > 0.0) scale_sq = d2;
        }
        effective = std::max(1e-30, scale_sq * 1e-22);
    }

    Mesh result;
    result.vertices = mesh.vertices;
    for (const auto& f : mesh.faces) {
        if (f[0] >= mesh.vertices.size() ||
            f[1] >= mesh.vertices.size() ||
            f[2] >= mesh.vertices.size())
            continue;  // skip malformed face
        const auto& v0 = mesh.vertices[f[0]];
        const auto& v1 = mesh.vertices[f[1]];
        const auto& v2 = mesh.vertices[f[2]];
        double area = internal::triangle_area(
            {v0[0], v0[1], v0[2]},
            {v1[0], v1[1], v1[2]},
            {v2[0], v2[1], v2[2]});
        if (area >= effective)
            result.faces.push_back(f);
    }
    return result;
}

Mesh remove_duplicate_faces(const Mesh& mesh) {
    struct FaceHash {
        size_t operator()(const std::array<uint32_t, 3>& f) const {
            size_t h = std::hash<uint32_t>{}(f[0]);
            h ^= std::hash<uint32_t>{}(f[1]) + 0x9e3779b9 + (h << 6) + (h >> 2);
            h ^= std::hash<uint32_t>{}(f[2]) + 0x9e3779b9 + (h << 6) + (h >> 2);
            return h;
        }
    };

    Mesh result;
    result.vertices = mesh.vertices;
    std::unordered_set<std::array<uint32_t, 3>, FaceHash> seen;
    for (const auto& f : mesh.faces) {
        std::array<uint32_t, 3> key = {f[0], f[1], f[2]};
        std::sort(key.begin(), key.end());
        if (seen.insert(key).second)
            result.faces.push_back(f);
    }
    return result;
}

DegenerateResult remove_degenerate(const Mesh& mesh, double area_threshold) {
    DegenerateResult dr{};

    Mesh after_degenerate = remove_degenerate_faces(mesh, area_threshold);
    dr.degenerate_removed = static_cast<uint32_t>(
        mesh.faces.size() - after_degenerate.faces.size());

    Mesh after_duplicate = remove_duplicate_faces(after_degenerate);
    dr.duplicate_removed = static_cast<uint32_t>(
        after_degenerate.faces.size() - after_duplicate.faces.size());

    // Compact vertices: keep only those referenced by remaining faces
    const auto& faces = after_duplicate.faces;
    const size_t old_vertex_count = after_duplicate.vertices.size();

    std::vector<uint32_t> remap(old_vertex_count, UINT32_MAX);
    for (const auto& f : faces) {
        if (f[0] >= remap.size() || f[1] >= remap.size() || f[2] >= remap.size())
            continue;  // skip malformed face during compaction
        remap[f[0]] = 0;
        remap[f[1]] = 0;
        remap[f[2]] = 0;
    }

    uint32_t new_index = 0;
    Mesh compacted;
    for (size_t i = 0; i < old_vertex_count; ++i) {
        if (remap[i] != UINT32_MAX) {
            remap[i] = new_index++;
            compacted.vertices.push_back(after_duplicate.vertices[i]);
        }
    }

    compacted.faces.reserve(faces.size());
    for (const auto& f : faces)
        compacted.faces.push_back({remap[f[0]], remap[f[1]], remap[f[2]]});

    dr.isolated_vertices_removed = static_cast<uint32_t>(
        old_vertex_count - compacted.vertices.size());
    dr.mesh = std::move(compacted);

    return dr;
}

} // namespace meshseal::stages
