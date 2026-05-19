#include <catch2/catch_test_macros.hpp>
#include <meshseal/meshseal.h>
#include "../meshseal/src/stages/holes.h"
#include "helpers.h"

using namespace meshseal;
// Topology helpers are provided by helpers.h (meshseal::test namespace).
using meshseal::test::boundary_edge_count;
using meshseal::test::is_consistently_wound;
using meshseal::test::max_edge_valence;
using meshseal::test::signed_volume;

TEST_CASE("Holes: open tetrahedron gets its base filled", "[holes]") {
    Mesh m = test::open_tetrahedron();
    REQUIRE(m.faces.size() == 3u);

    auto r = stages::fill_holes(m);
    CHECK(r.holes_filled == 1u);
    CHECK(r.faces_added >= 1u);
    CHECK(r.mesh.faces.size() > 3u);
}

TEST_CASE("Holes: open cube gets its rectangular hole filled", "[holes]") {
    Mesh m = test::open_cube();
    REQUIRE(m.faces.size() == 10u); // 12 - 2 top faces

    auto r = stages::fill_holes(m);
    CHECK(r.holes_filled == 1u);
    // Square hole (4-edge boundary loop) requires at least 2 triangles minimum,
    // and 4 triangles when using fan-from-centroid (which is used for loops ≤ max_fan_size).
    CHECK(r.faces_added >= 2u);
}

TEST_CASE("Holes: closed cube has no holes to fill", "[holes]") {
    Mesh m = test::unit_cube();

    auto r = stages::fill_holes(m);
    CHECK(r.holes_filled == 0u);
    CHECK(r.faces_added == 0u);
    CHECK(r.mesh.faces.size() == 12u);
}

TEST_CASE("Holes: empty mesh does not crash", "[holes]") {
    Mesh empty;
    stages::HoleFillResult r;
    REQUIRE_NOTHROW(r = stages::fill_holes(empty));
    CHECK(r.mesh.vertices.empty());
    CHECK(r.mesh.faces.empty());
    CHECK(r.holes_filled == 0u);
    CHECK(r.faces_added == 0u);
}

TEST_CASE("Holes: result has more faces than input after filling open tetrahedron", "[holes]") {
    Mesh m = test::open_tetrahedron();
    const auto faces_before = m.faces.size();

    auto r = stages::fill_holes(m);
    CHECK(r.mesh.faces.size() > faces_before);
}

TEST_CASE("Holes: filled open tetrahedron has no remaining boundary edges", "[holes]") {
    auto r = stages::fill_holes(test::open_tetrahedron());
    CHECK(boundary_edge_count(r.mesh) == 0u);
}

TEST_CASE("Holes: filled open cube has no remaining boundary edges", "[holes]") {
    auto r = stages::fill_holes(test::open_cube());
    CHECK(boundary_edge_count(r.mesh) == 0u);
}

TEST_CASE("Holes: mesh with two separate holes, both filled", "[holes]") {
    // Build a cube missing both the top (+Z) and bottom (-Z) faces → 2 holes.
    Mesh m = test::unit_cube();
    // Remove bottom (-Z) faces: indices 0,1 → {0,2,1},{0,3,2}
    // Remove top (+Z) faces: indices 2,3 → {4,5,6},{4,6,7}
    m.faces.erase(m.faces.begin(), m.faces.begin() + 4); // remove first 4 faces
    REQUIRE(m.faces.size() == 8u);
    REQUIRE(boundary_edge_count(m) == 8u); // 4 bottom + 4 top boundary edges

    auto r = stages::fill_holes(m);
    CHECK(r.holes_filled == 2u);
    CHECK(r.faces_added >= 4u); // at least 2 tris per hole
    CHECK(boundary_edge_count(r.mesh) == 0u);
}

TEST_CASE("Holes: filled open tetrahedron has consistent outward winding", "[holes]") {
    // boundary_edge_count==0 only proves closure; it cannot detect patch faces with
    // inward-facing normals.  This test independently verifies that fill triangles
    // are inserted with winding consistent with the surrounding mesh.
    auto r = stages::fill_holes(test::open_tetrahedron());
    CHECK(is_consistently_wound(r.mesh));
    CHECK(signed_volume(r.mesh) > 0.0);
    // Verify the patch region did not create a non-manifold edge
    CHECK(max_edge_valence(r.mesh) == 2u);
}

TEST_CASE("Holes: filled open cube has consistent outward winding", "[holes]") {
    auto r = stages::fill_holes(test::open_cube());
    CHECK(is_consistently_wound(r.mesh));
    CHECK(signed_volume(r.mesh) > 0.0);
    CHECK(max_edge_valence(r.mesh) == 2u);
}
