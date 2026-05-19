#include <catch2/catch_test_macros.hpp>
#include <meshseal/meshseal.h>
#include "../meshseal/src/stages/soup.h"
#include "helpers.h"

using namespace meshseal;
using meshseal::test::boundary_edge_count;
using meshseal::test::is_consistently_wound;
using meshseal::test::signed_volume;
using meshseal::test::max_edge_valence;

TEST_CASE("SoupReconstruct: triangle soup is recognised and reconstructed", "[soup_reconstruct]") {
    Mesh m = test::cube_triangle_soup();
    REQUIRE(m.vertices.size() == 36u);
    REQUIRE(m.faces.size() == 12u);

    auto r = stages::reconstruct_soup(m);
    CHECK(r.was_needed == true);
    CHECK(r.success == true);
    CHECK(!r.mesh.faces.empty());
    CHECK(r.faces_after == static_cast<uint32_t>(r.mesh.faces.size()));
    // The reconstructed cube must be a closed, correctly-wound manifold solid.
    CHECK(boundary_edge_count(r.mesh) == 0u);
    CHECK(is_consistently_wound(r.mesh));
    CHECK(signed_volume(r.mesh) > 0.0);
}

TEST_CASE("SoupReconstruct: already-manifold cube is a no-op", "[soup_reconstruct]") {
    Mesh m = test::unit_cube();
    const auto faces_before = m.faces.size();

    auto r = stages::reconstruct_soup(m);
    // A manifold single-component closed mesh does not need soup reconstruction.
    CHECK(r.was_needed == false);
    // When was_needed==false the stage skips reconstruction and returns the mesh unchanged.
    CHECK(r.mesh.faces.size() == faces_before);
    // Even though was_needed==false, the output must still be a valid closed manifold
    CHECK(boundary_edge_count(r.mesh) == 0u);
    CHECK(max_edge_valence(r.mesh) == 2u);
    CHECK(is_consistently_wound(r.mesh));
}

TEST_CASE("SoupReconstruct: empty mesh does not crash", "[soup_reconstruct]") {
    Mesh empty;
    stages::SoupResult r;
    REQUIRE_NOTHROW(r = stages::reconstruct_soup(empty));
    CHECK(r.mesh.vertices.empty());
    CHECK(r.mesh.faces.empty());
}

TEST_CASE("SoupReconstruct: faces_before matches input face count", "[soup_reconstruct]") {
    Mesh m = test::cube_triangle_soup();
    const auto input_faces = static_cast<uint32_t>(m.faces.size());

    auto r = stages::reconstruct_soup(m);
    CHECK(r.faces_before == input_faces);
}

TEST_CASE("SoupReconstruct: scrambled-normal soup is reconstructed to valid solid", "[soup_reconstruct]") {
    // All 12 face windings are flipped relative to cube_triangle_soup (normals point
    // inward).  The stage pre-processes with weld + orient before calling manifold,
    // so it must still produce a valid closed solid despite the scrambled normals.
    Mesh m = test::cube_triangle_soup();
    for (auto& f : m.faces) std::swap(f[1], f[2]); // flip every normal inward
    REQUIRE(m.vertices.size() == 36u);
    REQUIRE(m.faces.size() == 12u);

    auto r = stages::reconstruct_soup(m);
    CHECK(r.was_needed == true);
    CHECK(r.success == true);
    CHECK(!r.mesh.faces.empty());
    CHECK(boundary_edge_count(r.mesh) == 0u);
    CHECK(is_consistently_wound(r.mesh));
    CHECK(signed_volume(r.mesh) > 0.0);
}
