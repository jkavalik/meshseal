#include <catch2/catch_test_macros.hpp>
#include <meshseal/meshseal.h>
#include "../meshseal/src/stages/alpha_wrap.h"
#include "helpers.h"
#include <cmath>

using namespace meshseal;
using meshseal::test::boundary_edge_count;
using meshseal::test::signed_volume;
using meshseal::test::max_edge_valence;

TEST_CASE("alpha_wrap: closed cube wraps to a clean watertight solid",
          "[alpha_wrap]") {
    Mesh m = test::unit_cube();
    auto r = stages::alpha_wrap(m, /*alpha=*/-1.0, /*voxel_res=*/48);
    CHECK(r.success);
    REQUIRE(!r.mesh.faces.empty());
    CHECK(boundary_edge_count(r.mesh) == 0u);     // watertight
    CHECK(max_edge_valence(r.mesh) <= 2u);        // 2-manifold edges
    // The closing of a unit cube is the cube with its edges rounded by ~alpha
    // — volume stays in a sane band around 1.
    const double v = std::abs(signed_volume(r.mesh));
    CHECK(v > 0.3);
    CHECK(v < 4.0);
}

TEST_CASE("alpha_wrap: a gap narrower than 2*alpha is bridged",
          "[alpha_wrap]") {
    // open_cube has the +Z face removed — a 1x1 opening. With alpha = 0.6,
    // 2*alpha = 1.2 > 1, so the radius-alpha probe cannot fit through the
    // opening and the morphological closing seals it watertight.
    Mesh m = test::open_cube();
    REQUIRE(boundary_edge_count(m) > 0u);         // input is open
    auto r = stages::alpha_wrap(m, /*alpha=*/0.6, /*voxel_res=*/48);
    CHECK(r.success);
    REQUIRE(!r.mesh.faces.empty());
    CHECK(boundary_edge_count(r.mesh) == 0u);     // gap sealed -> watertight
    CHECK(max_edge_valence(r.mesh) <= 2u);
}

TEST_CASE("alpha_wrap: empty input fails cleanly", "[alpha_wrap]") {
    Mesh empty;
    auto r = stages::alpha_wrap(empty);
    CHECK_FALSE(r.success);
    CHECK(r.mesh.faces.empty());
}
