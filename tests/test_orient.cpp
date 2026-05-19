#include <catch2/catch_test_macros.hpp>
#include <meshseal/meshseal.h>
#include "../meshseal/src/stages/orient.h"
#include "helpers.h"
#include <utility>

using namespace meshseal;
// signed_volume and is_consistently_wound are provided by helpers.h (meshseal::test namespace).
// Bring them into file scope for use in CHECK expressions.
using meshseal::test::signed_volume;
using meshseal::test::is_consistently_wound;

TEST_CASE("Orient: all-inverted cube is fully corrected", "[orient]") {
    Mesh m = test::unit_cube_inverted();
    REQUIRE(m.faces.size() == 12u);

    auto r = stages::orient_mesh(m);
    CHECK(r.was_orientable == true);
    CHECK(r.components_flipped == 1u);
    CHECK(r.faces_flipped == 12u);
    CHECK(r.mesh.faces.size() == 12u);
    // directed-edge consistency: every edge appears once in each direction
    CHECK(is_consistently_wound(r.mesh));
    CHECK(signed_volume(r.mesh) > 0.0);
}

TEST_CASE("Orient: all-inverted tetrahedron is fully corrected", "[orient]") {
    Mesh m = test::unit_tetrahedron_inverted();
    REQUIRE(m.faces.size() == 4u);

    auto r = stages::orient_mesh(m);
    CHECK(r.was_orientable == true);
    CHECK(r.components_flipped == 1u);
    CHECK(r.faces_flipped == 4u);
    CHECK(r.mesh.faces.size() == 4u);
    CHECK(is_consistently_wound(r.mesh));
    CHECK(signed_volume(r.mesh) > 0.0);
}

TEST_CASE("Orient: already-correct cube has zero faces flipped", "[orient]") {
    Mesh m = test::unit_cube();
    REQUIRE(m.faces.size() == 12u);

    auto r = stages::orient_mesh(m);
    CHECK(r.faces_flipped == 0u);
    CHECK(r.mesh.faces.size() == 12u);
}

TEST_CASE("Orient: partially-flipped cube corrected", "[orient]") {
    Mesh m = test::cube_one_flipped();
    REQUIRE(m.faces.size() == 12u);

    auto r = stages::orient_mesh(m);
    CHECK(r.was_orientable == true);
    // BFS propagates consistent winding, then volume check flips the component if needed.
    // For a single-component mesh with 1 bad face, the algorithm may flip 1 face (correct
    // the outlier) or 11 faces (conform all others to the outlier then flip the whole
    // component). Both are correct outputs; verify the output mesh is correct rather than
    // counting flips.
    CHECK(r.mesh.faces.size() == 12u);
    CHECK((r.components_flipped == 0u || r.components_flipped == 1u));
    // Directed-edge consistency proves every individual face is correctly wound,
    // not just that the global aggregate volume happens to be positive.
    CHECK(is_consistently_wound(r.mesh));
    CHECK(signed_volume(r.mesh) > 0.0);
}

TEST_CASE("Orient: empty mesh does not crash", "[orient]") {
    Mesh empty;
    stages::OrientResult r;
    REQUIRE_NOTHROW(r = stages::orient_mesh(empty));
    CHECK(r.mesh.vertices.empty());
    CHECK(r.mesh.faces.empty());
    CHECK(r.faces_flipped == 0u);
}

TEST_CASE("Orient: two-component mesh with one inverted component", "[orient]") {
    // Two separate cubes: one correctly wound, one with all faces inverted.
    Mesh m = test::two_separate_cubes();
    REQUIRE(m.faces.size() == 24u);
    // Invert the second cube's faces (indices 12-23)
    for (size_t i = 12; i < m.faces.size(); ++i)
        std::swap(m.faces[i][1], m.faces[i][2]);

    auto r = stages::orient_mesh(m);
    CHECK(r.was_orientable == true);
    CHECK(r.mesh.faces.size() == 24u);
    CHECK(is_consistently_wound(r.mesh));
    CHECK(signed_volume(r.mesh) > 0.0);
    CHECK(r.components_flipped == 1u);
}

TEST_CASE("Orient: non-manifold edge sets was_orientable false", "[orient]") {
    // cube_with_fin has edge (0,1) shared by 3 faces — a genuine non-manifold edge.
    // orient_mesh must flag was_orientable=false, must not crash, and must not
    // remove any faces (orient only adjusts winding).
    Mesh m = test::cube_with_fin();
    REQUIRE(m.faces.size() == 13u);

    stages::OrientResult r;
    REQUIRE_NOTHROW(r = stages::orient_mesh(m));
    CHECK(r.was_orientable == false);
    CHECK(r.mesh.faces.size() == 13u);
    // orient must not partially modify a non-manifold mesh before bailing out
    CHECK(r.faces_flipped == 0u);
}
