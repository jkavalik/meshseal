#include "intersections.h"
#include <manifold/manifold.h>
#include <cstdint>
#include <exception>
#include <vector>

namespace meshseal::stages {

IntersectResult resolve_intersections(const Mesh& mesh) {
    IntersectResult result;
    result.faces_before = static_cast<uint32_t>(mesh.faces.size());
    result.faces_after  = result.faces_before;
    result.had_intersections = false;
    result.manifold_failed   = false;

    // Wrap the entire manifold-library interaction in a try/catch. Manifold
    // is documented to throw on some pathological inputs (very-large coords
    // and NaN cases that survive Phase 0); a throw here would propagate up
    // through every outer repair() frame and abort the whole pipeline.
    // Treat any throw the same as `manifold_failed`: return the original
    // mesh unchanged, caller's defect guards reject the no-op.
    try {

    // Build MeshGL from input Mesh
    manifold::MeshGL mesh_gl;
    mesh_gl.numProp = 3;

    mesh_gl.vertProperties.reserve(mesh.vertices.size() * 3);
    for (const auto& v : mesh.vertices) {
        mesh_gl.vertProperties.push_back(static_cast<float>(v[0]));
        mesh_gl.vertProperties.push_back(static_cast<float>(v[1]));
        mesh_gl.vertProperties.push_back(static_cast<float>(v[2]));
    }

    mesh_gl.triVerts.reserve(mesh.faces.size() * 3);
    for (const auto& f : mesh.faces) {
        mesh_gl.triVerts.push_back(static_cast<uint32_t>(f[0]));
        mesh_gl.triVerts.push_back(static_cast<uint32_t>(f[1]));
        mesh_gl.triVerts.push_back(static_cast<uint32_t>(f[2]));
    }

    // Construct Manifold
    manifold::Manifold m(mesh_gl);
    if (m.Status() != manifold::Manifold::Error::NoError) {
        result.manifold_failed = true;
        result.mesh = mesh;
        return result;
    }

    // Decompose into connected components, then union all of them.
    // A plain m + m self-union is the identity when m is already a valid
    // composite manifold (each point is already "inside" itself), so it
    // cannot remove intersection seams between distinct shells.
    // Computing the union of all decomposed components correctly resolves
    // interpenetrating shells.
    std::vector<manifold::Manifold> components = m.Decompose();
    manifold::Manifold repaired;
    if (components.size() <= 1) {
        // Single component: self-union is the only option.
        repaired = m + m;
    } else {
        repaired = components[0];
        for (size_t ci = 1; ci < components.size(); ++ci) {
            repaired = repaired + components[ci];
            if (repaired.Status() != manifold::Manifold::Error::NoError) {
                result.manifold_failed = true;
                result.mesh = mesh;
                return result;
            }
        }
    }

    if (repaired.Status() != manifold::Manifold::Error::NoError ||
        repaired.NumTri() == 0) {
        result.manifold_failed = true;
        result.mesh = mesh;
        return result;
    }

    // Extract result using double-precision output to avoid float32 truncation.
    manifold::MeshGL64 out = repaired.GetMeshGL64();

    if (out.triVerts.empty() || out.vertProperties.empty() || out.numProp < 3) {
        result.manifold_failed = true;
        result.mesh = mesh;
        return result;
    }

    // Convert MeshGL64 back to Mesh
    Mesh out_mesh;
    const uint64_t num_verts = static_cast<uint64_t>(out.vertProperties.size() / out.numProp);
    out_mesh.vertices.reserve(static_cast<size_t>(num_verts));
    for (uint64_t i = 0; i < num_verts; ++i) {
        uint64_t base = i * out.numProp;
        out_mesh.vertices.push_back({
            out.vertProperties[base + 0],
            out.vertProperties[base + 1],
            out.vertProperties[base + 2]
        });
    }

    const uint64_t num_tris = static_cast<uint64_t>(out.triVerts.size() / 3);
    out_mesh.faces.reserve(static_cast<size_t>(num_tris));
    for (uint64_t i = 0; i < num_tris; ++i) {
        uint64_t base = i * 3;
        out_mesh.faces.push_back({
            static_cast<uint32_t>(out.triVerts[base + 0]),
            static_cast<uint32_t>(out.triVerts[base + 1]),
            static_cast<uint32_t>(out.triVerts[base + 2])
        });
    }

    if (out_mesh.faces.empty() || out_mesh.vertices.empty()) {
        result.manifold_failed = true;
        result.mesh = mesh;
        return result;
    }

    result.faces_after       = static_cast<uint32_t>(out_mesh.faces.size());
    result.had_intersections = (result.faces_before != result.faces_after);
    result.mesh              = std::move(out_mesh);
    return result;
    } catch (const std::exception&) {
        result.manifold_failed = true;
        result.mesh = mesh;
        return result;
    } catch (...) {
        result.manifold_failed = true;
        result.mesh = mesh;
        return result;
    }
}

} // namespace meshseal::stages

