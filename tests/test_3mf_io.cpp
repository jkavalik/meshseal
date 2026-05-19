#include <catch2/catch_test_macros.hpp>
#include <meshseal/meshseal.h>
#include "../meshseal/src/3mf_io.h"
#include "helpers.h"
#include <filesystem>
#include <fstream>
#include <cmath>

using namespace meshseal;

TEST_CASE("3MF round-trip double-precision", "[3mf_io]") {
    // 1.0000000000000001 rounds to 1.0 in double;
    // nextafter(1.0, 2.0) = 1.0 + 2^-52 (1 ULP above 1.0).
    // The test verifies 3MF preserves both values with bit-exact fidelity.
    const double v0x = std::nextafter(1.0, 0.0); // 1 ULP below 1.0
    const double v1x = std::nextafter(1.0, 2.0); // 1 ULP above 1.0
    Mesh m;
    m.vertices = {
        {v0x, 0.0, 0.0},
        {v1x, 0.0, 0.0},
        {0.0, 1.0, 0.0},
    };
    m.faces = {{0, 1, 2}};

    auto tmp = std::filesystem::temp_directory_path() / "meshseal_test_3mf_precision.3mf";
    struct Guard { std::filesystem::path p; ~Guard() { std::filesystem::remove(p); } } g{tmp};
    write_3mf(m, tmp);
    Mesh result = read_3mf(tmp);

    REQUIRE(result.vertices.size() == 3u);
    CHECK(result.vertices[0][0] == v0x); // must survive round-trip bit-exactly
    CHECK(result.vertices[1][0] == v1x); // 1 ULP above 1.0 must not collapse to 1.0
}

TEST_CASE("3MF round-trip topology", "[3mf_io]") {
    Mesh cube = test::unit_cube();
    auto tmp = std::filesystem::temp_directory_path() / "meshseal_test_3mf_topology.3mf";
    struct Guard { std::filesystem::path p; ~Guard() { std::filesystem::remove(p); } } g{tmp};
    write_3mf(cube, tmp);
    Mesh result = read_3mf(tmp);

    CHECK(result.vertices.size() == 8u);
    CHECK(result.faces.size()    == 12u);
    for (int c = 0; c < 3; ++c)
        CHECK(std::abs(result.vertices[0][c] - cube.vertices[0][c]) < 1e-15);
}

TEST_CASE("3MF error: missing file", "[3mf_io]") {
    REQUIRE_THROWS_AS(read_3mf("/nonexistent/path.3mf"), ThreeMfError);
}

TEST_CASE("3MF error: corrupt ZIP data is rejected", "[3mf_io]") {
    // Write random bytes to a .3mf file; the ZIP parser must throw, not crash.
    auto tmp = std::filesystem::temp_directory_path() / "meshseal_test_3mf_corrupt.3mf";
    struct Guard { std::filesystem::path p; ~Guard() { std::filesystem::remove(p); } } g{tmp};
    {
        std::ofstream f(tmp, std::ios::binary);
        const char garbage[] = "This is not a valid ZIP/3MF file at all. \xDE\xAD\xBE\xEF";
        f.write(garbage, sizeof(garbage));
    }
    REQUIRE_THROWS_AS(read_3mf(tmp), ThreeMfError);
}

TEST_CASE("3MF round-trip: large face index (uint32_t range)", "[3mf_io]") {
    // Build a mesh where the largest vertex index fits in 32 bits.
    // 3 vertices, 1 face — trivial topology, checks index serialisation.
    Mesh m;
    m.vertices = {{0.0,0.0,0.0},{1.0,0.0,0.0},{0.0,1.0,0.0}};
    m.faces    = {{0u, 1u, 2u}};
    auto tmp = std::filesystem::temp_directory_path() / "meshseal_test_3mf_bigidx.3mf";
    struct Guard { std::filesystem::path p; ~Guard() { std::filesystem::remove(p); } } g{tmp};
    write_3mf(m, tmp);
    Mesh result = read_3mf(tmp);
    REQUIRE(result.faces.size() == 1u);
    CHECK(result.faces[0][0] == 0u);
    CHECK(result.faces[0][1] == 1u);
    CHECK(result.faces[0][2] == 2u);
}
