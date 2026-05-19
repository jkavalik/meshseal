#include <catch2/catch_test_macros.hpp>
#include <meshseal/meshseal.h>
#include "../meshseal/src/stages/nm_edge.h"
#include "helpers.h"

using namespace meshseal;
using meshseal::test::max_edge_valence;
using meshseal::test::boundary_edge_count;
using meshseal::test::is_consistently_wound;

TEST_CASE("NmEdge: fin face is removed from cube_with_fin", "[nm_edge]") {
    Mesh m = test::cube_with_fin();
    REQUIRE(m.faces.size() == 13u); // 12 cube + 1 fin

    auto r = stages::fix_non_manifold_edges(m);
    CHECK(r.had_non_manifold_edges == true);
    CHECK(r.faces_removed >= 1u);
    CHECK(r.edges_fixed >= 1u);
    // Exactly the fin must be gone; no cube face should be collateral damage.
    CHECK(r.mesh.faces.size() == 12u);
    // The result must be 2-manifold (every edge shared by exactly 2 faces).
    CHECK(max_edge_valence(r.mesh) == 2u);
    // Removing the fin must not open a boundary or break winding
    CHECK(boundary_edge_count(r.mesh) == 0u);
    CHECK(is_consistently_wound(r.mesh));
}

TEST_CASE("NmEdge: clean cube reports no non-manifold edges", "[nm_edge]") {
    Mesh m = test::unit_cube();

    auto r = stages::fix_non_manifold_edges(m);
    CHECK(r.had_non_manifold_edges == false);
    CHECK(r.faces_removed == 0u);
    CHECK(r.mesh.faces.size() == 12u);
}

TEST_CASE("NmEdge: empty mesh does not crash", "[nm_edge]") {
    Mesh empty;
    stages::NmEdgeResult r;
    REQUIRE_NOTHROW(r = stages::fix_non_manifold_edges(empty));
    CHECK(r.mesh.vertices.empty());
    CHECK(r.mesh.faces.empty());
    CHECK(r.had_non_manifold_edges == false);
    CHECK(r.faces_removed == 0u);
}

TEST_CASE("NmEdge: face count decreases after fixing fin", "[nm_edge]") {
    Mesh m = test::cube_with_fin();
    const auto faces_before = m.faces.size();

    auto r = stages::fix_non_manifold_edges(m);
    CHECK(r.mesh.faces.size() < faces_before);
}

TEST_CASE("NmEdge: multiple fins on separate edges both removed", "[nm_edge]") {
    // Build a cube with two independent fins on two different edges.
    // Edge (0,1) and edge (2,3) each get an extra fin triangle.
    Mesh m = test::unit_cube();
    // Fin 1: extra triangle on edge (0,1) — same as cube_with_fin helper
    m.vertices.push_back({0.5, -0.5, 0.0});  // index 8
    m.faces.push_back({0, 1, 8u});
    // Fin 2: extra triangle on edge (2,3) — cube has directed edges through those verts
    m.vertices.push_back({0.5, 1.5, 0.0});  // index 9
    m.faces.push_back({2, 3, 9u});
    REQUIRE(m.faces.size() == 14u);

    auto r = stages::fix_non_manifold_edges(m);
    CHECK(r.had_non_manifold_edges == true);
    CHECK(r.edges_fixed >= 2u);
    CHECK(r.faces_removed >= 2u);
    // Both fins gone; all 12 cube faces remain — not just ≤ 12 which would also
    // pass if the stage accidentally removed a cube face along with the fins.
    CHECK(r.mesh.faces.size() == 12u);
    CHECK(max_edge_valence(r.mesh) == 2u);
    CHECK(boundary_edge_count(r.mesh) == 0u);
    CHECK(is_consistently_wound(r.mesh));
}
