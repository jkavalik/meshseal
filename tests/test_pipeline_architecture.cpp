#include <catch2/catch_test_macros.hpp>
#include "../meshseal/include/meshseal/meshseal.h"
#include "helpers.h"
#include <cmath>

using namespace meshseal;
using namespace meshseal::test;

// ─── Test 1 ──────────────────────────────────────────────────────────────────
// SOUP routing: face_count ≤ 3 → silently discarded (not a failure)
// ─────────────────────────────────────────────────────────────────────────────
TEST_CASE("Pipeline: single-triangle soup component is discarded silently", "[pipeline][soup]")
{
    Mesh m;
    m.vertices = {{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}};
    m.faces    = {{0u, 1u, 2u}};

    // All stages enabled (default RepairOptions)
    RepairOptions opts;
    auto result = repair(m, opts);

    // A 1-face component has face_count ≤ 3 → SOUP → dropped; no volume content
    CHECK(result.mesh.faces.empty() == true);
    // Silent discard is not a pipeline failure
    CHECK(result.partial_failure == false);
    // Exactly one component outcome, classified as dropped_zero_volume with a
    // soup-noise rationale (1 face is below the soup-aggregation threshold).
    REQUIRE(result.component_outcomes.size() == 1u);
    bool found = false;
    for (const auto& co : result.component_outcomes) {
        if (co.status == ComponentStatus::dropped_zero_volume &&
            co.reason.find("soup") != std::string::npos) {
            found = true;
            break;
        }
    }
    CHECK(found);
}

// ─── Test 2 ──────────────────────────────────────────────────────────────────
// SOUP routing: open_ratio == 1.0 (all edges boundary) → soup_reconstruct
// ─────────────────────────────────────────────────────────────────────────────
TEST_CASE("Pipeline: cube_soup routes to soup_reconstruct and produces closed output", "[pipeline][soup]")
{
    // cube_triangle_soup() produces 12 independent triangles (cube faces
    // with no shared vertex indices).  Every edge is a boundary edge,
    // so open_ratio == 1.0, routing the component to the SOUP queue.
    Mesh m = cube_triangle_soup();

    RepairOptions opts;
    auto result = repair(m, opts);

    REQUIRE(result.mesh.faces.empty() == false);
    CHECK(boundary_edge_count(result.mesh) == 0u);
    CHECK(is_consistently_wound(result.mesh) == true);
    CHECK(signed_volume(result.mesh) > 0.5); // should recover a cube-like volume
}

// ─── Test 3 ──────────────────────────────────────────────────────────────────
// OPEN routing: flat coplanar patch — open_ratio = 0.5 (exactly at the planarity
// SOUP threshold, but NOT strictly greater-than), so it routes to OPEN.
// fill_holes closes the boundary loop, producing a zero-volume closed patch.
// analyze_shells max_abs_vol==0 guard skips the volume filter, so the closed
// zero-volume surface passes through — it is not a solid (is_volume=false).
// ─────────────────────────────────────────────────────────────────────────────
TEST_CASE("Pipeline: flat open triangle patch produces no solid output", "[pipeline][open]")
{
    // 4 coplanar triangles (all Z=0) forming a 2×2 grid.
    // open_ratio = 8 boundary / 16 total = 0.5 exactly.
    // Because the condition is open_ratio > 0.50 (strict), the patch routes to
    // OPEN, not SOUP.  fill_holes closes the border loop, producing a zero-volume
    // closed surface.  analyze_shells skips the volume filter when max_abs_vol==0,
    // so faces survive — but signed_volume==0 means is_volume==false (not a solid).
    Mesh m;
    m.vertices = {
        {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {2.0, 0.0, 0.0},
        {0.0, 1.0, 0.0}, {1.0, 1.0, 0.0}, {2.0, 1.0, 0.0},
        {0.0, 2.0, 0.0}, {1.0, 2.0, 0.0}, {2.0, 2.0, 0.0}
    };
    m.faces = {
        {0u, 1u, 3u}, {1u, 4u, 3u},
        {1u, 2u, 4u}, {2u, 5u, 4u},
        {3u, 4u, 6u}, {4u, 7u, 6u},
        {4u, 5u, 7u}, {5u, 8u, 7u}
    };

    RepairOptions opts;
    auto result = repair(m, opts);

    // The closed flat surface has no volumetric content: is_volume must be false
    // and the signed volume must be effectively zero.
    CHECK(result.is_volume == false);
    CHECK(std::abs(signed_volume(result.mesh)) < 1e-8);
}

// ─── Test 4 ──────────────────────────────────────────────────────────────────
// OPEN routing: fill_holes succeeds → component promoted to CLOSED queue
// ─────────────────────────────────────────────────────────────────────────────
TEST_CASE("Pipeline: open mesh gets hole filled and promoted to closed", "[pipeline][open]")
{
    // open_cube() = unit cube minus the two top (+Z) faces: 8 verts, 10 faces,
    // 1 open square boundary loop (4 boundary edges).
    // open_ratio > 0 → OPEN queue → fill_holes should close it → CLOSED queue.
    Mesh m = open_cube();

    RepairOptions opts;
    opts.weld          = false;
    opts.degenerate    = true;
    opts.orient        = true;
    opts.non_manifold  = true;
    opts.holes         = true;
    opts.shells        = true;
    opts.intersections = false;
    opts.thin_features = false;
    opts.soup_reconstruct = false;

    auto result = repair(m, opts);

    CHECK(boundary_edge_count(result.mesh) == 0u);
    CHECK(is_consistently_wound(result.mesh) == true);
    CHECK(result.watertight == true);
    CHECK(result.is_volume  == true);
    CHECK(result.mesh.faces.size() > 10u); // hole-fill added at least one face
}

// ─── Test 5 ──────────────────────────────────────────────────────────────────
// NO_BOUNDARY routing: already-closed mesh goes to CLOSED queue directly;
// no modifications expected on a geometrically valid input.
// ─────────────────────────────────────────────────────────────────────────────
TEST_CASE("Pipeline: valid closed mesh passes through without fill_holes", "[pipeline][no_boundary]")
{
    Mesh m = unit_cube(); // 8 verts, 12 faces, outward normals, no boundary edges
    REQUIRE(m.faces.size() == 12u);

    // All stages enabled — nothing to repair, pipeline should be a no-op
    RepairOptions opts;
    auto result = repair(m, opts);

    CHECK(result.watertight == true);
    CHECK(result.is_volume  == true);
    CHECK(boundary_edge_count(result.mesh) == 0u);
    CHECK(result.mesh.faces.size() == 12u); // no holes filled, no modifications
    CHECK(result.partial_failure == false);
}

// ─── Test 6 ──────────────────────────────────────────────────────────────────
// Phase 10R: partial_failure = true when UNREPAIRABLE PILE is non-empty.
// open_cube with holes=false: OPEN component cannot be closed → lands in pile.
// ─────────────────────────────────────────────────────────────────────────────
TEST_CASE("Pipeline: partial_failure is true when pile has components", "[pipeline][partial_failure]")
{
    Mesh m = open_cube();

    RepairOptions opts;
    opts.weld          = true;
    opts.degenerate    = true;
    opts.orient        = true;
    opts.non_manifold  = true;
    opts.holes         = false; // hole-filling disabled → open component cannot close
    opts.shells        = true;
    opts.intersections = false;
    opts.thin_features = false;
    opts.soup_reconstruct = false;

    auto result = repair(m, opts);

    CHECK(result.partial_failure == true);
    CHECK(result.mesh.faces.empty() == false); // pile faces are preserved in output
}

// ─── Test 7 ──────────────────────────────────────────────────────────────────
// Phase 9R: containment check — inner shell dropped with dropped_contained status
// ─────────────────────────────────────────────────────────────────────────────
TEST_CASE("Pipeline: contained inner shell is dropped with dropped_contained status", "[pipeline][shells]")
{
    // nested_cubes(): outer [0,4]³ (vol 64) + inner [1,3]³ (vol 8), both manifold.
    // The inner cube is fully contained by the outer; it should be dropped.
    Mesh m = nested_cubes();

    RepairOptions opts;
    opts.intersections = false;  // disable so manifold's boolean union does not
                                  // subsume the inner cube before analyze_shells
                                  // can classify it as contained
    auto result = repair(m, opts);

    CHECK(result.watertight == true);
    CHECK(result.is_volume  == true);
    CHECK(result.component_count == 1u);          // only outer shell survives
    CHECK(signed_volume(result.mesh) > 63.0);     // outer cube vol = 64
    CHECK(result.partial_failure == false);

    bool found_contained = false;
    for (const auto& co : result.component_outcomes) {
        if (co.status == ComponentStatus::dropped_contained) {
            found_contained = true;
            break;
        }
    }
    CHECK(found_contained);
}

// ─── Test 8 ──────────────────────────────────────────────────────────────────
// Phase 9R: 3 chain-overlapping shells merged via boolean union → single solid
// ─────────────────────────────────────────────────────────────────────────────
TEST_CASE("Pipeline: three chain-overlapping boxes merged into single solid", "[pipeline][intersections]")
{
    // three_chain_overlapping_boxes():
    //   A=[0,2]³, B=[1,3]×[0,2]×[0,2], C=[2,4]×[0,2]×[0,2]
    // Boolean union = [0,4]×[0,2]×[0,2], volume = 16.
    Mesh m = three_chain_overlapping_boxes();

    RepairOptions opts;
    auto result = repair(m, opts);

    CHECK(result.watertight == true);
    CHECK(result.is_volume  == true);
    CHECK(result.component_count == 1u);
    CHECK(signed_volume(result.mesh) > 15.5);
    CHECK(signed_volume(result.mesh) < 16.5);
    CHECK(boundary_edge_count(result.mesh) == 0u);
    CHECK(is_consistently_wound(result.mesh) == true);
}

// ─── Test 9 ──────────────────────────────────────────────────────────────────
// Phase 8R: thin fin excision on a closed component — fin removed, cube intact
// ─────────────────────────────────────────────────────────────────────────────
TEST_CASE("Pipeline: thin fin excision on closed component triggers fill_holes re-run", "[pipeline][thin_features]")
{
    Mesh m = unit_cube(); // 12 faces, closed manifold
    REQUIRE(m.faces.size() == 12u);

    // Append a zero-thickness fin: two back-to-back coplanar triangles using
    // existing cube vertices 0, 1, 2.  The pair is self-contained (their
    // directed edges pair up) so they add no new boundary edges of their own.
    m.faces.push_back({0u, 1u, 2u}); // fin face — forward winding
    m.faces.push_back({0u, 2u, 1u}); // fin face — reverse winding (zero thickness)
    // m now has 14 faces

    RepairOptions opts;
    opts.weld          = true;
    opts.degenerate    = true;
    opts.orient        = true;
    opts.non_manifold  = true;
    opts.holes         = true;
    opts.shells        = true;
    opts.intersections = false;
    opts.thin_features = true;
    opts.soup_reconstruct = false;

    auto result = repair(m, opts);

    CHECK(result.watertight == true);
    CHECK(boundary_edge_count(result.mesh) == 0u);
    CHECK(result.mesh.faces.size() == 12u); // fin removed; original 12 cube faces remain
}
