#include <catch2/catch_test_macros.hpp>
#include <meshseal/meshseal.h>
#include "../meshseal/src/stages/degenerate.h"
#include "helpers.h"

using namespace meshseal;

TEST_CASE("Degenerate: zero-area face is removed", "[degenerate]") {
    Mesh m = test::tetrahedron_with_degenerate();
    REQUIRE(m.faces.size() == 5u); // 4 valid + 1 degenerate

    auto r = stages::remove_degenerate(m);
    CHECK(r.degenerate_removed >= 1u);
    CHECK(r.mesh.faces.size() == 4u);
}

TEST_CASE("Degenerate: duplicate face is removed", "[degenerate]") {
    Mesh m = test::tetrahedron_with_duplicate();
    REQUIRE(m.faces.size() == 5u); // 4 valid + 1 duplicate

    auto r = stages::remove_degenerate(m);
    CHECK(r.duplicate_removed >= 1u);
    CHECK(r.mesh.faces.size() == 4u);
}

TEST_CASE("Degenerate: clean mesh has no faces removed", "[degenerate]") {
    Mesh m = test::unit_tetrahedron();
    REQUIRE(m.faces.size() == 4u);

    auto r = stages::remove_degenerate(m);
    CHECK(r.degenerate_removed == 0u);
    CHECK(r.duplicate_removed == 0u);
    CHECK(r.mesh.faces.size() == 4u);
}

TEST_CASE("Degenerate: empty mesh does not crash", "[degenerate]") {
    Mesh empty;
    stages::DegenerateResult r;
    REQUIRE_NOTHROW(r = stages::remove_degenerate(empty));
    CHECK(r.mesh.vertices.empty());
    CHECK(r.mesh.faces.empty());
    CHECK(r.degenerate_removed == 0u);
    CHECK(r.duplicate_removed == 0u);
}

TEST_CASE("Degenerate: isolated vertex is compacted", "[degenerate]") {
    Mesh m = test::unit_tetrahedron();
    // Add exactly one unreferenced vertex
    m.vertices.push_back({99.0, 99.0, 99.0});
    REQUIRE(m.vertices.size() == 5u);

    auto r = stages::remove_degenerate(m);
    // Exactly 1 isolated vertex was added; the count must be exactly 1, not just ≥ 1.
    CHECK(r.isolated_vertices_removed == 1u);
    CHECK(r.mesh.vertices.size() == 4u);
    // After compaction all face indices must still be in bounds.
    for (const auto& f : r.mesh.faces)
        for (int k = 0; k < 3; ++k)
            CHECK(f[k] < static_cast<uint32_t>(r.mesh.vertices.size()));
}

TEST_CASE("Degenerate: face with two identical indices is removed", "[degenerate]") {
    // {0, 0, 1} — two indices are the same; cross product is identically zero regardless
    // of where the vertices are. This pattern arises from edge-collapse operations.
    Mesh m = test::unit_tetrahedron();
    m.faces.push_back({0u, 0u, 1u}); // degenerate: v[0] == v[1]
    REQUIRE(m.faces.size() == 5u);

    auto r = stages::remove_degenerate(m);
    CHECK(r.degenerate_removed >= 1u);
    CHECK(r.mesh.faces.size() == 4u);
}

TEST_CASE("Degenerate: reversed duplicate is removed", "[degenerate]") {
    // {0,2,1} is already face 0 of the tetrahedron; {0,1,2} is its mirror (opposite winding).
    // A reversed duplicate is also a duplicate face and should be removed.
    Mesh m = test::unit_tetrahedron();
    const auto f0 = m.faces[0]; // {0,2,1}
    m.faces.push_back({f0[0], f0[2], f0[1]}); // reversed: {0,1,2}
    REQUIRE(m.faces.size() == 5u);

    auto r = stages::remove_degenerate(m);
    CHECK(r.duplicate_removed >= 1u);
    CHECK(r.mesh.faces.size() == 4u);
}

TEST_CASE("Degenerate: collinear face (distinct indices, zero area) is removed", "[degenerate]") {
    // Three distinct vertices collinear on the x-axis: their cross product is
    // identically zero, so area == 0.  This tests the area-based removal path
    // in remove_degenerate_faces (the index-equality path in weld would not catch it).
    Mesh m;
    m.vertices = {
        {0.0, 0.0, 0.0}, // 0
        {1.0, 0.0, 0.0}, // 1
        {2.0, 0.0, 0.0}, // 2 — collinear with 0 and 1
        {0.0, 0.0, 1.0}, // 3 — off-axis apex for a valid companion face
    };
    m.faces = {
        {0u, 1u, 2u}, // collinear: area = 0
        {0u, 1u, 3u}, // valid triangle
    };

    auto r = stages::remove_degenerate(m);
    CHECK(r.degenerate_removed == 1u);
    CHECK(r.mesh.faces.size() == 1u);
    // Vertex 2 is now unreferenced and must be compacted away.
    CHECK(r.isolated_vertices_removed == 1u);
    CHECK(r.mesh.vertices.size() == 3u);
    // All remaining face indices must be in bounds.
    for (const auto& f : r.mesh.faces)
        for (int k = 0; k < 3; ++k)
            CHECK(f[k] < static_cast<uint32_t>(r.mesh.vertices.size()));
}
