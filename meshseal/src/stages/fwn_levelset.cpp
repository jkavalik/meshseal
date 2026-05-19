#include "fwn_levelset.h"
#include "../internal/vec3.h"
#include "../internal/voxel_grid.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstring>
#include <limits>
#include <memory>
#include <string>
#include <unordered_map>
#include <manifold/common.h>
#include <manifold/manifold.h>

namespace meshseal::stages {

using internal::Vec3d;

namespace {

// Per-triangle solid angle contribution to the generalized winding number
// at query point p. From Jacobson, Kavan, Sorkine-Hornung (SIGGRAPH 2013),
// equation 6, using the robust atan2 form.
//
// Returns the SIGNED solid angle subtended by triangle (a,b,c) viewed from
// p, in steradians. Sign matches the triangle winding (right-hand rule).
//
// Sum over all triangles and divide by 4π to get the generalized winding
// number ω(p). For a closed orientable surface, ω is exactly +1 inside,
// 0 outside. For open/soup geometry, ω is fractional but the level set
// at 0.5 is a reasonable "best-guess inside" boundary.
inline double tri_solid_angle(const std::array<double, 3>& a,
                              const std::array<double, 3>& b,
                              const std::array<double, 3>& c,
                              const std::array<double, 3>& p)
{
    const double ax = a[0]-p[0], ay = a[1]-p[1], az = a[2]-p[2];
    const double bx = b[0]-p[0], by = b[1]-p[1], bz = b[2]-p[2];
    const double cx = c[0]-p[0], cy = c[1]-p[1], cz = c[2]-p[2];

    const double la = std::sqrt(ax*ax + ay*ay + az*az);
    const double lb = std::sqrt(bx*bx + by*by + bz*bz);
    const double lc = std::sqrt(cx*cx + cy*cy + cz*cz);

    // Degenerate: query point coincides with a vertex. Treat as 0.
    if (la < 1e-30 || lb < 1e-30 || lc < 1e-30) return 0.0;

    // Numerator: det([a-p, b-p, c-p]) = a·(b×c).
    const double det =
        ax * (by*cz - bz*cy) -
        ay * (bx*cz - bz*cx) +
        az * (bx*cy - by*cx);

    // Denominator: product term as per the robust atan2 form.
    const double dot_ab = ax*bx + ay*by + az*bz;
    const double dot_bc = bx*cx + by*cy + bz*cz;
    const double dot_ca = cx*ax + cy*ay + cz*az;
    const double denom = la*lb*lc + dot_ab*lc + dot_bc*la + dot_ca*lb;

    // Solid angle in steradians; 2·atan2 to get the full signed range.
    return 2.0 * std::atan2(det, denom);
}

} // anonymous namespace

FwnLevelSetResult fwn_levelset(const Mesh& mesh,
                               int    samples_per_axis,
                               double wn_threshold)
{
    FwnLevelSetResult result;
    result.faces_in = static_cast<uint32_t>(mesh.faces.size());

    if (mesh.faces.empty() || mesh.vertices.empty()) {
        result.reason = "empty input";
        return result;
    }

    // ---- bbox + padded bounds ----
    std::array<double, 3> lo = mesh.vertices[0];
    std::array<double, 3> hi = mesh.vertices[0];
    for (const auto& v : mesh.vertices) {
        for (int k = 0; k < 3; ++k) {
            lo[k] = std::min(lo[k], v[k]);
            hi[k] = std::max(hi[k], v[k]);
        }
    }
    const double ext_x = hi[0] - lo[0];
    const double ext_y = hi[1] - lo[1];
    const double ext_z = hi[2] - lo[2];
    const double diag = std::sqrt(ext_x*ext_x + ext_y*ext_y + ext_z*ext_z);
    if (diag <= 0.0) {
        result.reason = "degenerate bbox";
        return result;
    }
    // Pad bbox by 5 % so the isosurface doesn't clip the grid boundary.
    const double pad = 0.05 * diag;
    const manifold::vec3 bmin(lo[0]-pad, lo[1]-pad, lo[2]-pad);
    const manifold::vec3 bmax(hi[0]+pad, hi[1]+pad, hi[2]+pad);
    const manifold::Box bounds(bmin, bmax);

    const double max_ext = std::max({ext_x, ext_y, ext_z});
    if (samples_per_axis < 8) samples_per_axis = 8;
    const double edge_length = max_ext / static_cast<double>(samples_per_axis);

    // Pre-pack triangles into contiguous flat arrays to avoid indirection
    // in the inner loop of the SDF lambda (called once per marching-cubes
    // sample, typically 64³ = 262 K times).
    const size_t nF = mesh.faces.size();
    std::vector<std::array<double, 9>> tris(nF);
    for (size_t i = 0; i < nF; ++i) {
        const auto& f = mesh.faces[i];
        const auto& va = mesh.vertices[f[0]];
        const auto& vb = mesh.vertices[f[1]];
        const auto& vc = mesh.vertices[f[2]];
        tris[i] = {va[0], va[1], va[2],
                   vb[0], vb[1], vb[2],
                   vc[0], vc[1], vc[2]};
    }

    // SDF: s(p) = wn_threshold - ω(p).
    // Inside (ω > threshold) ⇒ s < 0 ; outside ⇒ s > 0 ; surface at s = 0.
    auto sdf = [tris, wn_threshold](manifold::vec3 mp) -> double {
        const std::array<double, 3> p = {mp.x, mp.y, mp.z};
        double sum = 0.0;
        for (const auto& t : tris) {
            sum += tri_solid_angle(
                {t[0], t[1], t[2]},
                {t[3], t[4], t[5]},
                {t[6], t[7], t[8]},
                p);
        }
        constexpr double kFourPi = 12.566370614359172;
        const double wn = sum / kFourPi;
        return wn_threshold - wn;
    };

    // Run marching cubes.
    manifold::Manifold m;
    try {
        m = manifold::Manifold::LevelSet(sdf, bounds, edge_length, /*level=*/0.0);
    } catch (const std::exception& e) {
        result.reason = std::string("LevelSet threw: ") + e.what();
        return result;
    }

    if (m.Status() != manifold::Manifold::Error::NoError || m.NumTri() == 0) {
        result.reason = "LevelSet produced empty mesh";
        return result;
    }

    // Convert back to Mesh.
    manifold::MeshGL out_gl = m.GetMeshGL();
    if (out_gl.triVerts.empty() || out_gl.vertProperties.empty() || out_gl.numProp < 3) {
        result.reason = "LevelSet produced invalid MeshGL";
        return result;
    }
    // MeshGL splits vertices per-property. Two triangles meeting at the
    // same 3D position can have different vertex indices if e.g. they
    // belong to different "groups". When we collapse to plain Mesh and
    // (later) cast to STL float32, those splits become spurious non-
    // manifold edges. Build a position-keyed weld so the returned Mesh
    // has one vertex per unique position.
    //
    // Use bit-exact float32 key (same encoding repair.cpp's float32-compat
    // uses): cast position component to float, reinterpret as uint32_t.
    // This catches all "would collide in STL" duplicates and ONLY those —
    // genuinely distinct positions stay distinct.
    struct F32Key {
        std::uint32_t x, y, z;
        bool operator==(const F32Key& o) const {
            return x == o.x && y == o.y && z == o.z;
        }
    };
    struct F32Hash {
        std::size_t operator()(const F32Key& k) const {
            std::size_t h = k.x;
            h ^= (std::size_t)k.y * 2654435761u + 0x9e3779b9u + (h << 6) + (h >> 2);
            h ^= (std::size_t)k.z * 2654435761u + 0x9e3779b9u + (h << 6) + (h >> 2);
            return h;
        }
    };
    auto pack = [](float f) -> std::uint32_t {
        std::uint32_t u;
        std::memcpy(&u, &f, sizeof(u));
        return u;
    };

    Mesh out;
    const uint32_t nv = static_cast<uint32_t>(out_gl.vertProperties.size() / out_gl.numProp);
    std::vector<uint32_t> remap(nv);
    std::unordered_map<F32Key, uint32_t, F32Hash> dedup;
    dedup.reserve(nv);
    for (uint32_t vi = 0; vi < nv; ++vi) {
        const uint32_t base = vi * out_gl.numProp;
        const float fx = out_gl.vertProperties[base + 0];
        const float fy = out_gl.vertProperties[base + 1];
        const float fz = out_gl.vertProperties[base + 2];
        F32Key k{pack(fx), pack(fy), pack(fz)};
        auto it = dedup.find(k);
        if (it == dedup.end()) {
            const uint32_t ni = static_cast<uint32_t>(out.vertices.size());
            out.vertices.push_back({
                static_cast<double>(fx),
                static_cast<double>(fy),
                static_cast<double>(fz)
            });
            dedup.emplace(k, ni);
            remap[vi] = ni;
        } else {
            remap[vi] = it->second;
        }
    }
    const uint32_t nt = static_cast<uint32_t>(out_gl.triVerts.size() / 3);
    out.faces.reserve(nt);
    for (uint32_t ti = 0; ti < nt; ++ti) {
        const uint32_t base = ti * 3;
        const uint32_t a = remap[out_gl.triVerts[base]];
        const uint32_t b = remap[out_gl.triVerts[base + 1]];
        const uint32_t c = remap[out_gl.triVerts[base + 2]];
        // Drop degenerate triangles created by the weld.
        if (a == b || b == c || c == a) continue;
        out.faces.push_back({a, b, c});
    }
    // Marching-cubes output can contain genuine NM edges at pinch points
    // between separately-extracted shells. Tried `fix_non_manifold_edges`
    // here to clean them — empirically regressed soup_seed sweep by 1
    // (27 → 28 failures), so we leave the FWN output unfiltered and rely
    // on the downstream Phase 8R + soup_reconstruct fallback to handle it.
    result.mesh = std::move(out);
    result.faces_out = static_cast<uint32_t>(result.mesh.faces.size());
    result.success = true;
    return result;
}

// ─── voxel_levelset ─────────────────────────────────────────────────────
// What `repair_algorithms_investigation.md` §1.3 actually prescribed for
// random-orientation triangle soups: voxelize input + flood-fill exterior +
// run marching cubes on the resulting occupancy field. Orientation-agnostic,
// which is what FWN is not.

namespace {

// Convert a manifold::Manifold to a Mesh with float32-bit-exact position
// weld + degenerate-face drop. Shared between fwn_levelset and
// voxel_levelset, since both produce dense MC output that needs the same
// post-treatment (MeshGL splits vertices per-property → spurious NM if not
// welded by position; float32 cast for STL preview must be exact).
Mesh manifold_to_welded_mesh(const manifold::Manifold& m, std::string& failure_reason) {
    Mesh out;
    if (m.Status() != manifold::Manifold::Error::NoError || m.NumTri() == 0) {
        failure_reason = "LevelSet produced empty mesh";
        return out;
    }
    manifold::MeshGL out_gl = m.GetMeshGL();
    if (out_gl.triVerts.empty() || out_gl.vertProperties.empty() || out_gl.numProp < 3) {
        failure_reason = "LevelSet produced invalid MeshGL";
        return out;
    }
    struct F32Key {
        std::uint32_t x, y, z;
        bool operator==(const F32Key& o) const {
            return x == o.x && y == o.y && z == o.z;
        }
    };
    struct F32Hash {
        std::size_t operator()(const F32Key& k) const {
            std::size_t h = k.x;
            h ^= (std::size_t)k.y * 2654435761u + 0x9e3779b9u + (h << 6) + (h >> 2);
            h ^= (std::size_t)k.z * 2654435761u + 0x9e3779b9u + (h << 6) + (h >> 2);
            return h;
        }
    };
    auto pack = [](float f) -> std::uint32_t {
        std::uint32_t u; std::memcpy(&u, &f, sizeof(u));
        return u;
    };
    const uint32_t nv = static_cast<uint32_t>(out_gl.vertProperties.size() / out_gl.numProp);
    std::vector<uint32_t> remap(nv);
    std::unordered_map<F32Key, uint32_t, F32Hash> dedup;
    dedup.reserve(nv);
    for (uint32_t vi = 0; vi < nv; ++vi) {
        const uint32_t base = vi * out_gl.numProp;
        const float fx = out_gl.vertProperties[base + 0];
        const float fy = out_gl.vertProperties[base + 1];
        const float fz = out_gl.vertProperties[base + 2];
        F32Key k{pack(fx), pack(fy), pack(fz)};
        auto it = dedup.find(k);
        if (it == dedup.end()) {
            const uint32_t ni = static_cast<uint32_t>(out.vertices.size());
            out.vertices.push_back({double(fx), double(fy), double(fz)});
            dedup.emplace(k, ni);
            remap[vi] = ni;
        } else {
            remap[vi] = it->second;
        }
    }
    const uint32_t nt = static_cast<uint32_t>(out_gl.triVerts.size() / 3);
    out.faces.reserve(nt);
    for (uint32_t ti = 0; ti < nt; ++ti) {
        const uint32_t a = remap[out_gl.triVerts[ti*3 + 0]];
        const uint32_t b = remap[out_gl.triVerts[ti*3 + 1]];
        const uint32_t c = remap[out_gl.triVerts[ti*3 + 2]];
        if (a == b || b == c || c == a) continue;
        out.faces.push_back({a, b, c});
    }
    return out;
}

} // anonymous namespace

FwnLevelSetResult voxel_levelset(const Mesh& mesh,
                                 int voxel_res,
                                 int samples_per_axis)
{
    FwnLevelSetResult result;
    result.faces_in = static_cast<uint32_t>(mesh.faces.size());

    if (mesh.faces.empty() || mesh.vertices.empty()) {
        result.reason = "empty input";
        return result;
    }

    // ---- bbox + padding ----
    std::array<double, 3> lo = mesh.vertices[0];
    std::array<double, 3> hi = mesh.vertices[0];
    for (const auto& v : mesh.vertices) {
        for (int k = 0; k < 3; ++k) {
            lo[k] = std::min(lo[k], v[k]);
            hi[k] = std::max(hi[k], v[k]);
        }
    }
    const double ext_x = hi[0] - lo[0];
    const double ext_y = hi[1] - lo[1];
    const double ext_z = hi[2] - lo[2];
    const double diag = std::sqrt(ext_x*ext_x + ext_y*ext_y + ext_z*ext_z);
    if (diag <= 0.0) {
        result.reason = "degenerate bbox";
        return result;
    }
    const double pad = 0.05 * diag;
    const manifold::vec3 bmin(lo[0]-pad, lo[1]-pad, lo[2]-pad);
    const manifold::vec3 bmax(hi[0]+pad, hi[1]+pad, hi[2]+pad);
    const manifold::Box bounds(bmin, bmax);
    const double max_ext = std::max({ext_x, ext_y, ext_z});
    if (samples_per_axis < 8) samples_per_axis = 8;
    if (voxel_res < 16)       voxel_res = 16;
    const double edge_length = max_ext / static_cast<double>(samples_per_axis);

    // ---- build voxel oracle ----
    auto grid = std::make_shared<internal::VoxelGrid>(
        lo, hi, static_cast<std::uint32_t>(voxel_res));
    for (const auto& f : mesh.faces) {
        grid->rasterize(mesh.vertices[f[0]], mesh.vertices[f[1]], mesh.vertices[f[2]]);
    }
    grid->flood_fill_exterior();

    // ---- early-out: is there any Inside cell at all? ----
    // If flood-fill marked everything Outside or Surface, there's no enclosed
    // pocket and LevelSet will return empty. Detect early and report cleanly
    // so the caller knows to keep the input as-is (the catastrophic-collapse
    // guard already does this for preserve_original components).
    {
        bool any_inside = false;
        const auto d = grid->dims();
        const double cs = grid->cell_size();
        // Sweep a small interior set of grid points to check (avoid scanning
        // every cell — for a 64³ grid that's 262K samples, expensive).
        // Sample on a 16³ subgrid.
        const std::uint32_t step = std::max<std::uint32_t>(1, d[0] / 16);
        for (std::uint32_t i = step; i < d[0] - step && !any_inside; i += step)
        for (std::uint32_t j = step; j < d[1] - step && !any_inside; j += step)
        for (std::uint32_t k = step; k < d[2] - step && !any_inside; k += step) {
            const std::array<double, 3> p = {
                lo[0] + (i + 0.5) * cs - 2.0 * cs,
                lo[1] + (j + 0.5) * cs - 2.0 * cs,
                lo[2] + (k + 0.5) * cs - 2.0 * cs
            };
            if (grid->sample(p) == internal::VoxelGrid::Label::Inside) any_inside = true;
        }
        if (!any_inside) {
            result.reason = "voxel oracle found no enclosed interior";
            return result;
        }
    }

    // ---- SDF: voxel-occupancy → ±1 with Surface as 0 ----
    // Inside  → -1.0 (s<0, inside the solid)
    // Outside → +1.0 (s>0, outside the solid)
    // Surface →  0.0 (on the boundary; level=0 picks it)
    // Marching cubes interpolates between corner SDF values, so a binary
    // ±1 SDF gives the cell-boundary-style "voxel cube" output, which is
    // exactly what voxelization + MC is supposed to produce. The Surface
    // = 0 case nudges the iso to sit on the marked band when triangulating.
    auto sdf = [grid](manifold::vec3 mp) -> double {
        const std::array<double, 3> p = {mp.x, mp.y, mp.z};
        using L = internal::VoxelGrid::Label;
        const L lab = grid->sample(p);
        if (lab == L::Inside)  return -1.0;
        if (lab == L::Outside) return  1.0;
        return 0.0; // Surface
    };

    manifold::Manifold m;
    try {
        m = manifold::Manifold::LevelSet(sdf, bounds, edge_length, /*level=*/0.0);
    } catch (const std::exception& e) {
        result.reason = std::string("LevelSet threw: ") + e.what();
        return result;
    }

    std::string fail;
    Mesh out = manifold_to_welded_mesh(m, fail);
    if (out.faces.empty()) {
        result.reason = fail;
        return result;
    }

    result.mesh = std::move(out);
    result.faces_out = static_cast<uint32_t>(result.mesh.faces.size());
    result.success = true;
    return result;
}

} // namespace meshseal::stages
