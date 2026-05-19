#include <catch2/catch_test_macros.hpp>
#include <meshseal/meshseal.h>
#include "../meshseal/src/stages/nm_vertex.h"
#include "helpers.h"
#include <algorithm>
#include <map>
#include <set>
#include <vector>

using namespace meshseal;
using meshseal::test::boundary_edge_count;
using meshseal::test::is_consistently_wound;

// NOTE: This helper defines ring-adjacency as "two faces incident to vertex v share
// any third vertex != v". This is correct for simple fixtures like bowtie_tetrahedra.
// However it can produce false negatives for vertices where two separate fans
// accidentally share a third vertex — use only on the fixtures defined in this file.
// Helper: count vertices whose one-ring is not a single connected disk.
// Returns the number of non-manifold vertices.
static uint32_t count_non_manifold_vertices(const Mesh& m) {
    // For each vertex v, collect its incident faces and build an adjacency graph
    // over those faces: two faces in the ring of v are adjacent iff they share
    // a neighbour of v (i.e. they share the ring edge through v).
    std::map<uint32_t, std::vector<uint32_t>> v_to_faces;
    for (uint32_t fi = 0; fi < static_cast<uint32_t>(m.faces.size()); ++fi)
        for (int k = 0; k < 3; ++k)
            v_to_faces[m.faces[fi][k]].push_back(fi);

    uint32_t nm = 0;
    for (const auto& [v, fids] : v_to_faces) {
        if (fids.empty()) continue;

        // Build adjacency: fi and fj are ring-adjacent iff they share a vertex != v.
        std::map<uint32_t, std::vector<uint32_t>> adj;
        for (uint32_t fi : fids) {
            for (uint32_t fj : fids) {
                if (fi >= fj) continue;
                std::set<uint32_t> si(m.faces[fi].begin(), m.faces[fi].end());
                std::set<uint32_t> sj(m.faces[fj].begin(), m.faces[fj].end());
                for (uint32_t x : si) {
                    if (x != v && sj.count(x)) {
                        adj[fi].push_back(fj);
                        adj[fj].push_back(fi);
                        break;
                    }
                }
            }
        }

        // BFS: if the one-ring is not connected → non-manifold vertex.
        std::set<uint32_t> visited;
        std::vector<uint32_t> queue = {fids[0]};
        visited.insert(fids[0]);
        while (!queue.empty()) {
            uint32_t cur = queue.back(); queue.pop_back();
            for (uint32_t nb : adj[cur])
                if (!visited.count(nb)) { visited.insert(nb); queue.push_back(nb); }
        }
        if (visited.size() < fids.size()) ++nm;
    }
    return nm;
}

TEST_CASE("NmVertex: bowtie vertex is split", "[nm_vertex]") {
    Mesh m = test::bowtie_tetrahedra();
    const auto verts_before = m.vertices.size();

    auto r = stages::split_non_manifold_vertices(m);
    CHECK(r.vertices_split >= 1u);
    CHECK(r.vertices_added >= 1u);
    CHECK(r.mesh.vertices.size() > verts_before);
}

TEST_CASE("NmVertex: clean cube has no vertices split", "[nm_vertex]") {
    Mesh m = test::unit_cube();

    auto r = stages::split_non_manifold_vertices(m);
    CHECK(r.vertices_split == 0u);
    CHECK(r.vertices_added == 0u);
    CHECK(r.mesh.vertices.size() == m.vertices.size());
}

TEST_CASE("NmVertex: empty mesh does not crash", "[nm_vertex]") {
    Mesh empty;
    stages::NmVertexResult r;
    REQUIRE_NOTHROW(r = stages::split_non_manifold_vertices(empty));
    CHECK(r.mesh.vertices.empty());
    CHECK(r.mesh.faces.empty());
    CHECK(r.vertices_split == 0u);
    CHECK(r.vertices_added == 0u);
}

TEST_CASE("NmVertex: result has strictly more vertices than input for bowtie", "[nm_vertex]") {
    Mesh m = test::bowtie_tetrahedra();
    const auto verts_before = m.vertices.size();
    const auto faces_before = m.faces.size();

    auto r = stages::split_non_manifold_vertices(m);
    CHECK(r.mesh.vertices.size() > verts_before);
    // Face count must be preserved — only vertices are added, not faces
    CHECK(r.mesh.faces.size() == faces_before);
}

TEST_CASE("NmVertex: bowtie output has no remaining non-manifold vertices", "[nm_vertex]") {
    auto r = stages::split_non_manifold_vertices(test::bowtie_tetrahedra());
    CHECK(count_non_manifold_vertices(r.mesh) == 0u);
    // All face indices must be valid after the split
    for (const auto& f : r.mesh.faces)
        for (int k = 0; k < 3; ++k)
            CHECK(f[k] < static_cast<uint32_t>(r.mesh.vertices.size()));
    // Splitting a vertex must not introduce boundary edges or corrupt winding
    CHECK(boundary_edge_count(r.mesh) == 0u);
    CHECK(is_consistently_wound(r.mesh));
}
