#include <catch2/catch_test_macros.hpp>
#include <meshseal/meshseal.h>
#include "../meshseal/src/stages/shells.h"
#include "helpers.h"

using namespace meshseal;
using meshseal::test::boundary_edge_count;
using meshseal::test::is_consistently_wound;
using meshseal::test::max_edge_valence;

TEST_CASE("Shells: two equal-volume cubes both kept", "[shells]") {
    Mesh m = test::two_separate_cubes();

    auto r = stages::analyze_shells(m);
    CHECK(r.shells_kept == 2u);
    CHECK(r.shells_dropped == 0u);
    CHECK(r.shells.size() == 2u);
}

TEST_CASE("Shells: tiny cube dropped, big cube kept", "[shells]") {
    Mesh m = test::big_and_tiny_cube();

    auto r = stages::analyze_shells(m);
    CHECK(r.shells_dropped == 1u);
    CHECK(r.shells_kept == 1u);
}

TEST_CASE("Shells: single cube produces one shell", "[shells]") {
    Mesh m = test::unit_cube();

    auto r = stages::analyze_shells(m);
    CHECK(r.shells_kept == 1u);
    CHECK(r.shells_dropped == 0u);
    CHECK(r.shells.size() == 1u);
}

TEST_CASE("Shells: empty mesh does not crash", "[shells]") {
    Mesh empty;
    stages::ShellResult r;
    REQUIRE_NOTHROW(r = stages::analyze_shells(empty));
    CHECK(r.mesh.vertices.empty());
    CHECK(r.mesh.faces.empty());
    CHECK(r.shells_kept == 0u);
    CHECK(r.shells_dropped == 0u);
}

TEST_CASE("Shells: dropped shell is not in output mesh", "[shells]") {
    // After analyze_shells, the output mesh must contain only the kept shells.
    auto r = stages::analyze_shells(test::big_and_tiny_cube());
    // Big cube kept (≈volume 1), tiny cube dropped (volume ≈1e-9 < 1% of 1).
    CHECK(r.shells_kept == 1u);
    // The output mesh must have fewer faces than the 24-face combined input.
    CHECK(r.mesh.faces.size() < test::big_and_tiny_cube().faces.size());
    // The output mesh must be non-empty (big cube remains).
    CHECK(!r.mesh.faces.empty());
    // Verify the surviving big cube is still a valid closed manifold
    CHECK(boundary_edge_count(r.mesh) == 0u);
    CHECK(max_edge_valence(r.mesh) == 2u);
    CHECK(is_consistently_wound(r.mesh));
}

TEST_CASE("Shells: ShellInfo volume_threshold boundary - shell above 1% is kept", "[shells]") {
    // Build a big cube (volume 1) and a small cube with volume > 0.01 (above 1% threshold).
    // With vol_threshold=0.01, the condition to drop is |vol| < 0.01 * max_vol.
    // s=0.22 → vol = 0.22^3 = 0.010648 > 0.01 → must be kept.
    Mesh m = test::unit_cube(); // vol = 1.0
    const double s = 0.22; // s^3 = 0.010648 > 0.01
    Mesh small;
    small.vertices = {
        {5.0,   5.0,   5.0},   {5.0+s, 5.0,   5.0},
        {5.0+s, 5.0+s, 5.0},   {5.0,   5.0+s, 5.0},
        {5.0,   5.0,   5.0+s}, {5.0+s, 5.0,   5.0+s},
        {5.0+s, 5.0+s, 5.0+s}, {5.0,   5.0+s, 5.0+s}
    };
    small.faces = {
        {0,2,1},{0,3,2}, {4,5,6},{4,6,7},
        {0,1,5},{0,5,4}, {2,3,7},{2,7,6},
        {0,4,7},{0,7,3}, {1,2,6},{1,6,5}
    };
    const uint32_t off = static_cast<uint32_t>(m.vertices.size());
    for (auto& f : small.faces) f = {f[0]+off, f[1]+off, f[2]+off};
    m.vertices.insert(m.vertices.end(), small.vertices.begin(), small.vertices.end());
    m.faces.insert(m.faces.end(), small.faces.begin(), small.faces.end());

    auto r = stages::analyze_shells(m, 0.01);
    CHECK(r.shells_kept == 2u);
    CHECK(r.shells_dropped == 0u);
}
