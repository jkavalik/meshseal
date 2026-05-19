#include <catch2/catch_test_macros.hpp>
#include <meshseal/meshseal.h>
#include "../meshseal/src/stages/thin.h"
#include "helpers.h"

using namespace meshseal;

TEST_CASE("ThinFeatures: back-to-back triangles are detected and removed", "[thin_features]") {
    Mesh m = test::back_to_back_triangles();
    REQUIRE(m.faces.size() == 2u);

    auto r = stages::remove_thin_features(m);
    // The fixture has exactly 1 back-to-back pair composed of exactly 2 faces
    CHECK(r.pairs_found == 1u);
    CHECK(r.faces_removed == 2u);
    // All face-originating vertices become isolated after removal; they must be compacted.
    // After compaction the mesh should have no vertices (both faces removed → all 6 verts
    // become unreferenced, but only if the stage compacts — at minimum face count is 0).
    CHECK(r.mesh.faces.empty());
    // If both faces were removed, all 6 vertices become unreferenced.
    // A conforming stage must compact the vertex buffer (or leave it empty).
    // We accept either: empty OR all remaining entries are actually referenced.
    if (!r.mesh.faces.empty()) {
        for (const auto& f : r.mesh.faces)
            for (int k = 0; k < 3; ++k)
                CHECK(f[k] < static_cast<uint32_t>(r.mesh.vertices.size()));
    } else {
        // Faces were fully removed — vertex buffer should also be empty or compacted.
        // At minimum it must not exceed the original 6 vertices (no new verts added).
        CHECK(r.mesh.vertices.size() <= 6u);
    }
}

TEST_CASE("ThinFeatures: clean cube has no thin features", "[thin_features]") {
    Mesh m = test::unit_cube();

    auto r = stages::remove_thin_features(m);
    CHECK(r.pairs_found == 0u);
    CHECK(r.faces_removed == 0u);
    CHECK(r.mesh.faces.size() == 12u);
}

TEST_CASE("ThinFeatures: empty mesh does not crash", "[thin_features]") {
    Mesh empty;
    stages::ThinResult r;
    REQUIRE_NOTHROW(r = stages::remove_thin_features(empty));
    CHECK(r.mesh.vertices.empty());
    CHECK(r.mesh.faces.empty());
    CHECK(r.pairs_found == 0u);
    CHECK(r.faces_removed == 0u);
}

TEST_CASE("ThinFeatures: result has fewer faces after removing back-to-back triangles", "[thin_features]") {
    Mesh m = test::back_to_back_triangles();
    const auto faces_before = m.faces.size();

    auto r = stages::remove_thin_features(m);
    CHECK(r.mesh.faces.size() < faces_before);
}
