#pragma once
#include "../../include/meshseal/meshseal.h"
#include <cstdint>
#include <cstring>
#include <string>
#include <unordered_map>
#include <vector>
#include <manifold/manifold.h>

namespace meshseal::internal {

// Convert a manifold::Manifold to a Mesh with float32-bit-exact position
// weld + degenerate-face drop. Shared by every volumetric reconstructor
// (fwn_levelset, voxel_levelset, alpha_wrap): they all produce dense
// marching-cubes output that needs the same post-treatment — MeshGL splits
// vertices per-property, which yields spurious non-manifold edges unless the
// output is welded by exact float32 position, and the float32 cast is what
// an STL preview will see anyway.
inline Mesh manifold_to_welded_mesh(const manifold::Manifold& m,
                                    std::string& failure_reason) {
    Mesh out;
    if (m.Status() != manifold::Manifold::Error::NoError || m.NumTri() == 0) {
        failure_reason = "LevelSet produced empty mesh";
        return out;
    }
    manifold::MeshGL out_gl = m.GetMeshGL();
    if (out_gl.triVerts.empty() || out_gl.vertProperties.empty() ||
        out_gl.numProp < 3) {
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
    const std::uint32_t nv =
        static_cast<std::uint32_t>(out_gl.vertProperties.size() / out_gl.numProp);
    std::vector<std::uint32_t> remap(nv);
    std::unordered_map<F32Key, std::uint32_t, F32Hash> dedup;
    dedup.reserve(nv);
    for (std::uint32_t vi = 0; vi < nv; ++vi) {
        const std::uint32_t base = vi * out_gl.numProp;
        const float fx = out_gl.vertProperties[base + 0];
        const float fy = out_gl.vertProperties[base + 1];
        const float fz = out_gl.vertProperties[base + 2];
        F32Key k{pack(fx), pack(fy), pack(fz)};
        auto it = dedup.find(k);
        if (it == dedup.end()) {
            const std::uint32_t ni = static_cast<std::uint32_t>(out.vertices.size());
            out.vertices.push_back({double(fx), double(fy), double(fz)});
            dedup.emplace(k, ni);
            remap[vi] = ni;
        } else {
            remap[vi] = it->second;
        }
    }
    const std::uint32_t nt =
        static_cast<std::uint32_t>(out_gl.triVerts.size() / 3);
    out.faces.reserve(nt);
    for (std::uint32_t ti = 0; ti < nt; ++ti) {
        const std::uint32_t a = remap[out_gl.triVerts[ti*3 + 0]];
        const std::uint32_t b = remap[out_gl.triVerts[ti*3 + 1]];
        const std::uint32_t c = remap[out_gl.triVerts[ti*3 + 2]];
        if (a == b || b == c || c == a) continue;
        out.faces.push_back({a, b, c});
    }
    return out;
}

} // namespace meshseal::internal
