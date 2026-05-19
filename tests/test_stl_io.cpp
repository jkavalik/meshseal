#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <meshseal/meshseal.h>
#include "../meshseal/src/stl_io.h"
#include "helpers.h"
#include <filesystem>
#include <cstdio>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <vector>

using namespace meshseal;

TEST_CASE("STL round-trip preserves topology", "[stl_io]") {
    Mesh cube = test::unit_cube();
    auto bytes = write_stl_bytes(cube);
    Mesh result = read_stl_bytes(bytes.data(), bytes.size());

    // Binary STL stores 3 unshared vertices per triangle
    REQUIRE(result.vertices.size() == 36u);
    REQUIRE(result.faces.size() == 12u);
    for (int c = 0; c < 3; ++c) {
        CHECK(std::isfinite(result.vertices.front()[c]));
        CHECK(std::isfinite(result.vertices.back()[c]));
    }
}

TEST_CASE("STL round-trip file", "[stl_io]") {
    auto tmp = std::filesystem::temp_directory_path() / "meshseal_test_stl_roundtrip.stl";
    struct Guard { std::filesystem::path p; ~Guard() { std::filesystem::remove(p); } } g{tmp};
    write_stl(test::unit_cube(), tmp);
    Mesh result = read_stl(tmp);

    CHECK(result.faces.size() == 12u);
    for (auto& v : result.vertices) {
        for (int c = 0; c < 3; ++c)
            CHECK(std::isfinite(v[c]));
    }
}

TEST_CASE("STL validation: file too small", "[stl_io]") {
    uint8_t buf[10] = {};
    REQUIRE_THROWS_AS(read_stl_bytes(buf, sizeof(buf)), StlError);
}

TEST_CASE("STL validation: size mismatch", "[stl_io]") {
    // Header (80) + tri_count=1 (4) = 84 bytes, but zero triangle data.
    // Expected total = 84 + 1*50 = 134; actual = 84 → StlError.
    std::vector<uint8_t> buf(84, 0);
    buf[80] = 1; // tri_count = 1 (little-endian)
    REQUIRE_THROWS_AS(read_stl_bytes(buf.data(), buf.size()), StlError);
}

TEST_CASE("STL validation: NaN vertex", "[stl_io]") {
    Mesh tri;
    tri.vertices = {{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}};
    tri.faces    = {{0, 1, 2}};
    auto bytes = write_stl_bytes(tri);
    REQUIRE(bytes.size() == 134u);

    // Offset 96 = 80 (header) + 4 (tri_count) + 12 (normal) = first vertex x.
    // Write 0x7FC00000 (quiet NaN) in little-endian.
    constexpr uint32_t nan_bits = 0x7FC00000u;
    bytes[96] = static_cast<uint8_t>( nan_bits        & 0xFF);
    bytes[97] = static_cast<uint8_t>((nan_bits >>  8) & 0xFF);
    bytes[98] = static_cast<uint8_t>((nan_bits >> 16) & 0xFF);
    bytes[99] = static_cast<uint8_t>((nan_bits >> 24) & 0xFF);

    REQUIRE_THROWS_AS(read_stl_bytes(bytes.data(), bytes.size()), StlError);
}

TEST_CASE("STL valid-mesh pass-through: tetrahedron", "[stl_io]") {
    auto bytes = write_stl_bytes(test::unit_tetrahedron());
    Mesh result = read_stl_bytes(bytes.data(), bytes.size());
    CHECK(result.faces.size() == 4u);
    for (auto& v : result.vertices)
        for (int c = 0; c < 3; ++c)
            CHECK(std::isfinite(v[c]));
}

TEST_CASE("STL validation: Inf vertex is rejected", "[stl_io]") {
    Mesh tri;
    tri.vertices = {{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}};
    tri.faces    = {{0, 1, 2}};
    auto bytes = write_stl_bytes(tri);
    REQUIRE(bytes.size() == 134u);

    // Write +Inf (0x7F800000) into first vertex x (same offset as NaN test).
    constexpr uint32_t inf_bits = 0x7F800000u;
    bytes[96] = static_cast<uint8_t>( inf_bits        & 0xFF);
    bytes[97] = static_cast<uint8_t>((inf_bits >>  8) & 0xFF);
    bytes[98] = static_cast<uint8_t>((inf_bits >> 16) & 0xFF);
    bytes[99] = static_cast<uint8_t>((inf_bits >> 24) & 0xFF);

    REQUIRE_THROWS_AS(read_stl_bytes(bytes.data(), bytes.size()), StlError);
}

TEST_CASE("STL round-trip: vertex coordinates survive at float32 precision", "[stl_io]") {
    // Binary STL stores float32; verify that a value exactly representable in float32
    // round-trips without error (no extra precision loss beyond the format's inherent limit).
    Mesh tri;
    // Use vertices with exact float32 representations.
    tri.vertices = {{0.5f, 0.25f, 0.125f}, {1.0f, 0.0f, 0.0f}, {0.0f, 1.0f, 0.0f}};
    tri.faces    = {{0, 1, 2}};
    auto bytes  = write_stl_bytes(tri);
    Mesh result = read_stl_bytes(bytes.data(), bytes.size());
    REQUIRE(result.vertices.size() == 3u);
    CHECK(result.vertices[0][0] == static_cast<double>(0.5f));
    CHECK(result.vertices[0][1] == static_cast<double>(0.25f));
    CHECK(result.vertices[0][2] == static_cast<double>(0.125f));
}

TEST_CASE("STL round-trip: empty mesh writes and reads back as empty", "[stl_io]") {
    Mesh empty;
    auto bytes  = write_stl_bytes(empty);
    // Header (80) + tri_count=0 (4) = 84 bytes.
    REQUIRE(bytes.size() == 84u);
    Mesh result = read_stl_bytes(bytes.data(), bytes.size());
    CHECK(result.vertices.empty());
    CHECK(result.faces.empty());
}

TEST_CASE("STL read: ASCII STL is detected and parsed", "[stl_io]") {
    // A minimal ASCII STL with two triangles (a quad). ASCII STLs are common
    // in real-world inputs; the reader must not mistake them for binary.
    const char* ascii =
        "solid quad\n"
        "  facet normal 0 0 1\n"
        "    outer loop\n"
        "      vertex 0 0 0\n"
        "      vertex 1 0 0\n"
        "      vertex 1 1 0\n"
        "    endloop\n"
        "  endfacet\n"
        "  facet normal 0 0 1\n"
        "    outer loop\n"
        "      vertex 0 0 0\n"
        "      vertex 1 1 0\n"
        "      vertex 0 1 0\n"
        "    endloop\n"
        "  endfacet\n"
        "endsolid quad\n";
    const auto* data = reinterpret_cast<const uint8_t*>(ascii);
    Mesh m = read_stl_bytes(data, std::strlen(ascii));
    REQUIRE(m.faces.size() == 2u);
    REQUIRE(m.vertices.size() == 6u);
    CHECK(m.vertices[1][0] == 1.0);
    CHECK(m.vertices[5][1] == 1.0);
}

TEST_CASE("STL read: ASCII STL with scientific-notation coords", "[stl_io]") {
    // Real exporters emit scientific notation (e.g. 8.603526e-01).
    const char* ascii =
        "solid s\n"
        "facet normal -5.077787e-01 -4.420782e-02 8.603526e-01\n"
        "outer loop\n"
        "vertex 8.750794e-02 -3.810792e-07 4.295658e+01\n"
        "vertex 8.618489e-02 1.519706e-02 4.295658e+01\n"
        "vertex 0.0 0.0 0.0\n"
        "endloop\n"
        "endfacet\n"
        "endsolid s\n";
    const auto* data = reinterpret_cast<const uint8_t*>(ascii);
    Mesh m = read_stl_bytes(data, std::strlen(ascii));
    REQUIRE(m.faces.size() == 1u);
    CHECK(m.vertices[0][2] == Catch::Approx(42.95658));
    CHECK(m.vertices[1][1] == Catch::Approx(0.01519706));
}
