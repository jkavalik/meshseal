#include <catch2/catch_test_macros.hpp>
#include <meshseal/meshseal.h>
#include "../meshseal/src/stages/intersections.h"
#include "helpers.h"

using namespace meshseal;
using meshseal::test::boundary_edge_count;
using meshseal::test::is_consistently_wound;
using meshseal::test::max_edge_valence;
using meshseal::test::signed_volume;

TEST_CASE("Intersections: overlapping tetrahedra resolved without manifold failure", "[intersections]") {
    Mesh m = test::two_overlapping_tetrahedra();
    REQUIRE(m.faces.size() == 8u);

    stages::IntersectResult r;
    REQUIRE_NOTHROW(r = stages::resolve_intersections(m));
    CHECK(r.faces_before == 8u);
    CHECK(r.manifold_failed == false);
    // Result must be non-empty (self-union should yield a valid solid)
    CHECK(!r.mesh.faces.empty());
    CHECK(r.faces_after == static_cast<uint32_t>(r.mesh.faces.size()));
    // The boolean union of two overlapping convex solids must be a closed,
    // consistently-wound manifold with positive volume.
    CHECK(boundary_edge_count(r.mesh) == 0u);
    CHECK(is_consistently_wound(r.mesh));
    CHECK(signed_volume(r.mesh) > 0.0);
}

TEST_CASE("Intersections: clean cube is a no-op (no intersections detected)", "[intersections]") {
    Mesh m = test::unit_cube();
    REQUIRE(m.faces.size() == 12u);

    auto r = stages::resolve_intersections(m);
    CHECK(r.manifold_failed == false);
    // Clean single-component closed mesh: stage should not detect intersections
    CHECK(r.had_intersections == false);
    // Face count should be unchanged when no intersection resolution was needed
    CHECK(r.faces_before == 12u);
    CHECK(r.faces_after == 12u);
}

TEST_CASE("Intersections: empty mesh does not crash", "[intersections]") {
    Mesh empty;
    stages::IntersectResult r;
    REQUIRE_NOTHROW(r = stages::resolve_intersections(empty));
    CHECK(r.mesh.vertices.empty());
    CHECK(r.mesh.faces.empty());
    // Empty mesh may cause manifold to fail gracefully — that is acceptable
}

TEST_CASE("Intersections: faces_before and faces_after are consistent", "[intersections]") {
    Mesh m = test::unit_cube();
    const auto input_faces = static_cast<uint32_t>(m.faces.size());

    auto r = stages::resolve_intersections(m);
    CHECK(r.faces_before == input_faces);
    CHECK(r.faces_after == static_cast<uint32_t>(r.mesh.faces.size()));
}

// ─── Failure-mode regression tests ───────────────────────────────────────────
// These tests document the two categories of failure revealed by the benchmark:
//
//  • Open/non-manifold input  → manifold construction fails (NotManifold) →
//    stage correctly passes the mesh through unchanged. Responsible for the 13
//    real-world models (airplane_builtin, bird_bath, black_vase, lobster, …)
//    that still have self-intersections after repair: the stage cannot help them
//    because they arrive here with open boundary edges or non-manifold edge
//    topology that manifold's constructor rejects.
//
//  • Multi-component overlapping shells → stage works correctly when called,
//    but the pipeline guard `orig_has_open_boundaries` prevents it from running
//    for originally-open combined-defect meshes (kitchen_sink,
//    overlapping_shells_with_hole, clown, doorman, human, o3d_monkey …).
//    This is the pipeline-level bug; see the final test case below.

TEST_CASE("Intersections: open boundary mesh is returned unchanged (manifold rejects it)", "[intersections]") {
    // Root cause of the 13 NotManifold real-world failures.
    // An open mesh (boundary edges present) cannot be loaded into manifold; the
    // stage must preserve the original and set manifold_failed=true.
    Mesh m = test::open_cube();  // unit cube with top (+Z) face removed: 10 faces
    REQUIRE(m.faces.size() == 10u);

    auto r = stages::resolve_intersections(m);

    CHECK(r.manifold_failed == true);
    // Mesh returned bit-for-bit unchanged
    CHECK(r.mesh.faces.size()    == m.faces.size());
    CHECK(r.mesh.vertices.size() == m.vertices.size());
    CHECK(r.had_intersections    == false);
    CHECK(r.faces_before         == 10u);
    CHECK(r.faces_after          == 10u);
}

TEST_CASE("Intersections: non-manifold edge mesh is returned unchanged (manifold rejects it)", "[intersections]") {
    // Root cause of disjoint_fans and similar cases where nm_edge repair leaves
    // a mesh with edges shared by 3+ faces — manifold's constructor returns
    // NotManifold and the stage must pass the mesh through.
    Mesh m = test::cube_with_fin();  // unit cube + 1 fin face making one edge 3-valent
    REQUIRE(m.faces.size() == 13u);
    REQUIRE(max_edge_valence(m) == 3u);  // confirm the NM edge is present

    auto r = stages::resolve_intersections(m);

    CHECK(r.manifold_failed == true);
    CHECK(r.mesh.faces.size()    == m.faces.size());
    CHECK(r.mesh.vertices.size() == m.vertices.size());
    CHECK(r.had_intersections    == false);
    CHECK(r.faces_before         == 13u);
    CHECK(r.faces_after          == 13u);
}

TEST_CASE("Intersections: two overlapping closed cubes are merged by boolean union", "[intersections]") {
    // The stage itself works correctly for closed multi-component overlapping input.
    // This is the topology reached by meshes like kitchen_sink or
    // overlapping_shells_with_hole AFTER their holes have been filled — but the
    // pipeline guard currently prevents the stage from running for those cases
    // (see the pipeline regression test below).
    //
    // A=[0,1]^3, B=[0.5,1.5]×[0,1]×[0,1].  Overlap vol=0.5, union vol=1.5.
    Mesh m = test::two_overlapping_cubes();
    REQUIRE(m.faces.size() == 24u);

    auto r = stages::resolve_intersections(m);

    REQUIRE(r.manifold_failed == false);
    CHECK(boundary_edge_count(r.mesh) == 0u);
    CHECK(is_consistently_wound(r.mesh));
    CHECK(signed_volume(r.mesh) > 0.0);
    // Union of two unit cubes with 0.5 overlap = volume 1.5
    CHECK(signed_volume(r.mesh) > 1.45);
    CHECK(signed_volume(r.mesh) < 1.55);
    // Face count changes because the overlap region introduces new boundary edges
    CHECK(r.had_intersections == true);
}

TEST_CASE("Intersections: two non-overlapping shells leave had_intersections false", "[intersections]") {
    // When components do not intersect geometrically, manifold's union preserves
    // the total face count exactly.  had_intersections must remain false to avoid
    // false positives in stage notes.
    Mesh m = test::two_separate_cubes();  // one at origin, one at (10,10,10)
    REQUIRE(m.faces.size() == 24u);

    auto r = stages::resolve_intersections(m);

    REQUIRE(r.manifold_failed == false);
    CHECK(boundary_edge_count(r.mesh) == 0u);
    CHECK(is_consistently_wound(r.mesh));
    // Non-overlapping union: signed volume = sum of both unit cubes = 2.0
    CHECK(signed_volume(r.mesh) > 1.95);
    CHECK(signed_volume(r.mesh) < 2.05);
    // No intersection seams → face count must not change
    CHECK(r.had_intersections == false);
}

TEST_CASE("Intersections: nested shells — inner cube absorbed into outer", "[intersections]") {
    // Boolean union of two nested solids equals the outer solid, because every
    // point of the inner is already inside the outer.  Represents the
    // matryoshka_nested / nested_shells fixture patterns.
    //
    // Outer=[0,4]^3 (vol 64), Inner=[1,3]^3 (vol 8). Union vol = 64.
    Mesh m = test::nested_cubes();
    REQUIRE(m.faces.size() == 24u);

    auto r = stages::resolve_intersections(m);

    REQUIRE(r.manifold_failed == false);
    CHECK(boundary_edge_count(r.mesh) == 0u);
    CHECK(is_consistently_wound(r.mesh));
    // Union equals the outer cube: volume must be close to 64
    CHECK(signed_volume(r.mesh) > 63.5);
    CHECK(signed_volume(r.mesh) < 64.5);
    // Inner faces (12) are absorbed → face count changes
    CHECK(r.had_intersections == true);
}

TEST_CASE("Intersections: three chain-overlapping boxes merged into single solid", "[intersections]") {
    // Tests the Decompose() + iterative union loop for n_components > 2.
    // Represents real-world patterns with 3-6 overlapping components:
    // clown (4 comp), doorman (6 comp), human (4 comp), o3d_monkey (3 comp).
    //
    // A=[0,2]×[0,2]×[0,2], B=[1,3]×[0,2]×[0,2], C=[2,4]×[0,2]×[0,2].
    // A∩B=vol 4, B∩C=vol 4, A∩C=0. Union=[0,4]×[0,2]×[0,2]=vol 16.
    Mesh m = test::three_chain_overlapping_boxes();
    REQUIRE(m.faces.size() == 36u);  // 3 × 12

    auto r = stages::resolve_intersections(m);

    REQUIRE(r.manifold_failed == false);
    CHECK(boundary_edge_count(r.mesh) == 0u);
    CHECK(is_consistently_wound(r.mesh));
    // Volume must be the union (16), not the sum (24)
    CHECK(signed_volume(r.mesh) > 15.5);
    CHECK(signed_volume(r.mesh) < 16.5);
    CHECK(r.had_intersections == true);
}

// ─── Pipeline regression test ─────────────────────────────────────────────────

TEST_CASE("Intersections: pipeline merges overlapping shells after hole-filling (regression)", "[intersections]") {
    // BUG: The full repair() pipeline skips the intersection stage whenever the
    // *original* input had open boundary edges (`orig_has_open_boundaries=true`),
    // even after those boundaries have been closed by the holes stage.
    //
    // Fixture: open_cube [0,1]^3 (top face missing) + complete cube
    // [0.5,1.5]×[0,1]×[0,1].  After hole-filling, both shells are closed and
    // their overlap region [0.5,1]×[0,1]×[0,1] (vol 0.5) should be resolved by
    // the intersection stage.  The expected result is a single merged solid with
    // volume 1.5.
    //
    // Currently FAILS: the pipeline returns component_count=2 (two separate
    // watertight shells that still overlap in space) instead of 1 (merged solid),
    // because the guard `!orig_has_open_boundaries` in repair.cpp prevents the
    // intersection stage from running.
    Mesh m = test::overlapping_cubes_one_holed();
    REQUIRE(m.faces.size() == 22u);  // 10 (open cube) + 12 (closed cube)

    RepairOptions opts;
    opts.weld = opts.degenerate = opts.orient = opts.non_manifold =
        opts.holes = opts.shells = opts.intersections = true;

    auto result = repair(m, opts);

    CHECK(result.watertight == true);
    CHECK(result.is_volume  == true);
    // After merging the two shells the result must be a single solid, not two
    // separate overlapping components.  This CHECK currently fails.
    CHECK(result.component_count == 1);
    // Union volume 1.5 = 1 (cube A) + 1 (cube B) - 0.5 (overlap)
    CHECK(signed_volume(result.mesh) > 1.45);
    CHECK(signed_volume(result.mesh) < 1.55);
}
