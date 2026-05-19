#include <catch2/catch_test_macros.hpp>
#include <meshseal/meshseal.h>
#include "../meshseal/src/stages/weld.h"
#include "helpers.h"

using namespace meshseal;
using meshseal::test::is_consistently_wound;
using meshseal::test::max_edge_valence;
using meshseal::test::boundary_edge_count;

TEST_CASE("Weld: STL soup -> shared topology", "[weld]") {
    // Build an unshared vertex soup: 36 vertices, 12 faces (3 unique verts per tri)
    Mesh soup;
    Mesh cube = test::unit_cube();
    for (auto& face : cube.faces) {
        uint32_t base = static_cast<uint32_t>(soup.vertices.size());
        soup.vertices.push_back(cube.vertices[face[0]]);
        soup.vertices.push_back(cube.vertices[face[1]]);
        soup.vertices.push_back(cube.vertices[face[2]]);
        soup.faces.push_back({base, base + 1u, base + 2u});
    }
    REQUIRE(soup.vertices.size() == 36u);
    REQUIRE(soup.faces.size()    == 12u);

    auto r = stages::weld(soup, 1e-6);
    CHECK(r.vertices_before == 36u);
    CHECK(r.vertices_after  == 8u);
    CHECK(r.mesh.faces.size() == 12u);
    // After welding, the cube must be a closed 2-manifold with outward normals.
    CHECK(max_edge_valence(r.mesh) == 2u);
    CHECK(is_consistently_wound(r.mesh));
    CHECK(boundary_edge_count(r.mesh) == 0u);
}

TEST_CASE("Weld: degenerate faces dropped after weld", "[weld]") {
    Mesh m;
    m.vertices = {
        {0.0, 0.0, 0.0},  // 0
        {0.0, 0.0, 0.0},  // 1 — same position as 0
        {1.0, 0.0, 0.0},  // 2
        {0.0, 1.0, 0.0},  // 3
    };
    // Face 0 uses v0 and v1 (identical positions) → becomes degenerate after weld
    // Face 1 is a valid triangle
    m.faces = {{0, 1, 2}, {0, 2, 3}};

    auto r = stages::weld(m, 1e-6);
    // v0 and v1 merge; face {0,1,2} collapses to a line → dropped
    CHECK(r.mesh.faces.size() == 1u);
}

TEST_CASE("Weld: empty mesh", "[weld]") {
    Mesh empty;
    stages::WeldResult r;
    REQUIRE_NOTHROW(r = stages::weld(empty));
    CHECK(r.mesh.vertices.empty());
    CHECK(r.mesh.faces.empty());
}

TEST_CASE("Weld: already-shared mesh unchanged", "[weld]") {
    Mesh tet = test::unit_tetrahedron();
    auto r = stages::weld(tet, 1e-6);
    // No vertices should be merged — the tetrahedron is already welded.
    CHECK(r.vertices_after == r.vertices_before);
    CHECK(r.mesh.faces.size() == tet.faces.size());
}

TEST_CASE("Weld: tolerance edge case", "[weld]") {
    const double tol = 0.1;

    // Distance 0.05 < tol: the two close vertices should be merged,
    // collapsing the single face to degenerate (it gets dropped).
    Mesh m1;
    m1.vertices = {
        {0.00, 0.0, 0.0},
        {0.05, 0.0, 0.0},  // 0.05 < tol
        {0.00, 1.0, 0.0},
    };
    m1.faces = {{0, 1, 2}};
    auto r1 = stages::weld(m1, tol);
    CHECK(r1.vertices_after < 3u); // v0 and v1 merged

    // Distance 0.2 > tol: the two vertices should NOT be merged.
    Mesh m2;
    m2.vertices = {
        {0.00, 0.0, 0.0},
        {0.20, 0.0, 0.0},  // 0.2 > tol
        {0.00, 1.0, 0.0},
    };
    m2.faces = {{0, 1, 2}};
    auto r2 = stages::weld(m2, tol);
    CHECK(r2.vertices_after == 3u); // no merging

    // Distance exactly == tol: should be merged (comparison is <=)
    Mesh m3;
    m3.vertices = {
        {0.00, 0.0, 0.0},
        {0.10, 0.0, 0.0},  // distance == tol exactly
        {0.00, 1.0, 0.0},
    };
    m3.faces = {{0, 1, 2}};
    auto r3 = stages::weld(m3, tol);
    CHECK(r3.vertices_after < 3u); // v0 and v1 must merge at boundary
}

TEST_CASE("Weld: auto-tolerance (nullopt) welds cube soup correctly", "[weld]") {
    // nullopt triggers the scale-adaptive tolerance computation path:
    //   tolerance = max(1e-9, mean_edge_length * 1e-4, bbox_diagonal * 1e-7)
    // For a unit cube soup: diagonal ≈ 1.73, mean edge ≈ 1.0,
    // so tolerance ≈ max(1e-9, 1e-4, 1.73e-7) = 1e-4.  All 36 unshared vertices
    // must collapse to the 8 unique cube corners.
    Mesh soup = test::cube_triangle_soup();
    REQUIRE(soup.vertices.size() == 36u);

    auto r = stages::weld(soup, std::nullopt);
    CHECK(r.vertices_after == 8u);
    CHECK(r.mesh.faces.size() == 12u);
    CHECK(r.tolerance_used > 0.0); // confirms the auto path ran
    CHECK(max_edge_valence(r.mesh) == 2u);
    CHECK(is_consistently_wound(r.mesh));
    CHECK(boundary_edge_count(r.mesh) == 0u);
}
