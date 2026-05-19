#include <catch2/catch_test_macros.hpp>
#include <meshseal/meshseal.h>
#include "helpers.h"
#include <algorithm>

using namespace meshseal;
using meshseal::test::is_consistently_wound;
using meshseal::test::signed_volume;
using meshseal::test::boundary_edge_count;

TEST_CASE("repair: valid mesh passes through", "[e2e]") {
    auto result = repair(test::unit_cube());
    CHECK(result.watertight == true);
    CHECK(result.is_volume  == true);
    // Independent geometric checks (not pipeline self-assessment):
    CHECK(is_consistently_wound(result.mesh));
    CHECK(signed_volume(result.mesh) > 0.0);
}

TEST_CASE("repair: empty mesh", "[e2e]") {
    Mesh empty;
    RepairResult result;
    REQUIRE_NOTHROW(result = repair(empty));
    // Empty mesh has no edges → diagnostics considers it vacuously watertight
    CHECK(result.watertight == true);
}

TEST_CASE("repair: stages_applied contains weld", "[e2e]") {
    auto result = repair(test::unit_cube());
    bool has_weld = std::find(result.stages_applied.begin(),
                              result.stages_applied.end(),
                              "weld") != result.stages_applied.end();
    CHECK(has_weld);
}

TEST_CASE("repair: stage_times_ms populated", "[e2e]") {
    auto result = repair(test::unit_cube());
    CHECK(result.stage_times_ms.count("weld") > 0);
    CHECK(result.stage_times_ms.at("weld") >= 0.0);
}

TEST_CASE("repair: degenerate faces removed", "[e2e][degenerate]") {
    using namespace meshseal::test;
    auto mesh = tetrahedron_with_degenerate();
    auto result = meshseal::repair(mesh);
    // The degenerate face {0,0,0} is caught by weld (identical indices) or degenerate stage;
    // either way the output must have fewer faces than the input.
    // tetrahedron_with_degenerate has exactly 1 extra degenerate face — verify exactly 1 removed
    CHECK(result.mesh.faces.size() == mesh.faces.size() - 1);
}

TEST_CASE("repair: duplicate faces removed", "[e2e][degenerate]") {
    using namespace meshseal::test;
    auto mesh = tetrahedron_with_duplicate();
    auto result = meshseal::repair(mesh);
    // tetrahedron_with_duplicate has exactly 1 extra face
    CHECK(result.mesh.faces.size() == mesh.faces.size() - 1);
    bool has_duplicate_note = false;
    for (const auto& n : result.notes)
        if (n.find("duplicate") != std::string::npos) has_duplicate_note = true;
    CHECK(has_duplicate_note);
}

TEST_CASE("repair: inverted normals fixed", "[e2e][orient]") {
    auto mesh = meshseal::test::unit_cube_inverted();
    auto result = meshseal::repair(mesh);
    CHECK(result.watertight);
    CHECK(result.is_volume);
    bool has_orient_note = false;
    for (const auto& n : result.notes)
        if (n.find("orient") != std::string::npos || n.find("flipped") != std::string::npos)
            has_orient_note = true;
    CHECK(has_orient_note);
    // Pipeline self-assessment is insufficient: result.watertight and result.is_volume
    // are computed by the pipeline's own diagnostic code, so a bug in both the repair
    // stage and the diagnostic would still pass those checks. Verify independently.
    CHECK(is_consistently_wound(result.mesh));
    CHECK(signed_volume(result.mesh) > 0.0);
}

TEST_CASE("repair: already-correct normals not flipped", "[e2e][orient]") {
    auto result = meshseal::repair(meshseal::test::unit_cube());
    CHECK(result.watertight);
    CHECK(result.is_volume);
    // If winding is already correct, faces_flipped should be 0.
    // Verify independently — pipeline self-assessment (watertight/is_volume) is not enough.
    CHECK(is_consistently_wound(result.mesh));
    CHECK(signed_volume(result.mesh) > 0.0);
}

TEST_CASE("repair: bowtie vertex is split", "[e2e][nm_vertex]") {
    auto mesh = meshseal::test::bowtie_tetrahedra();
    // Disable stages that re-merge coincident vertices (manifold) or may drop one of the two tets
    meshseal::RepairOptions opts;
    opts.holes            = false;
    opts.shells           = false;
    opts.intersections    = false;
    opts.soup_reconstruct = false;
    auto result = meshseal::repair(mesh, opts);
    // After repair: more vertices than before (vertex 0 was split)
    CHECK(result.mesh.vertices.size() > mesh.vertices.size());
    // nm_vertex should be in stages_applied
    bool has_nm_stage = false;
    for (const auto& s : result.stages_applied)
        if (s == "nm_vertex") has_nm_stage = true;
    CHECK(has_nm_stage);
    // After the full repair pipeline the mesh must be a valid closed solid
    CHECK(boundary_edge_count(result.mesh) == 0u);
    CHECK(meshseal::test::max_edge_valence(result.mesh) == 2u);
    CHECK(is_consistently_wound(result.mesh));
}

TEST_CASE("repair: manifold mesh unchanged by nm_vertex", "[e2e][nm_vertex]") {
    auto mesh = meshseal::test::unit_cube();
    auto result = meshseal::repair(mesh);
    // Cube has no non-manifold vertices, vertex count must stay the same
    CHECK(result.mesh.vertices.size() == mesh.vertices.size());
}

TEST_CASE("repair: fin face removed", "[e2e][nm_edge]") {
    auto mesh = meshseal::test::cube_with_fin();
    auto result = meshseal::repair(mesh);
    // After repair, the extra fin face should be gone (strictly fewer faces)
    CHECK(result.mesh.faces.size() < mesh.faces.size());
    // The result should be watertight (pipeline self-assessment)
    CHECK(result.watertight);
    // Independent topology verification — pipeline diagnostics are not sufficient
    CHECK(boundary_edge_count(result.mesh) == 0u);
    CHECK(meshseal::test::max_edge_valence(result.mesh) == 2u);
    CHECK(is_consistently_wound(result.mesh));
}

TEST_CASE("repair: manifold mesh unchanged by nm_edge", "[e2e][nm_edge]") {
    auto result = meshseal::repair(meshseal::test::unit_cube());
    CHECK(result.watertight);
    CHECK(result.mesh.faces.size() == meshseal::test::unit_cube().faces.size());
}

TEST_CASE("repair: open cube filled", "[e2e][holes]") {
    auto result = meshseal::repair(meshseal::test::open_cube());
    CHECK(result.watertight);
    bool note_found = false;
    for (const auto& n : result.notes)
        if (n.find("holes") != std::string::npos) note_found = true;
    CHECK(note_found);
    // Verify closure independently — no boundary edges must remain.
    CHECK(boundary_edge_count(result.mesh) == 0u);
}

TEST_CASE("repair: open tetrahedron filled", "[e2e][holes]") {
    auto result = meshseal::repair(meshseal::test::open_tetrahedron());
    CHECK(result.watertight);
    CHECK(result.is_volume);
    CHECK(boundary_edge_count(result.mesh) == 0u);
    CHECK(is_consistently_wound(result.mesh));
    CHECK(signed_volume(result.mesh) > 0.0);
}

TEST_CASE("repair: closed cube not modified by holes", "[e2e][holes]") {
    auto mesh = meshseal::test::unit_cube();
    auto result = meshseal::repair(mesh);
    CHECK(result.mesh.faces.size() == mesh.faces.size());
}

TEST_CASE("repair: two separate cubes kept", "[e2e][shells]") {
    auto result = meshseal::repair(meshseal::test::two_separate_cubes());
    // Both shells have equal volume so neither is dropped
    CHECK(result.mesh.faces.size() == meshseal::test::two_separate_cubes().faces.size());
    CHECK(result.component_count == 2);
}

TEST_CASE("repair: tiny shell dropped", "[e2e][shells]") {
    auto mesh = meshseal::test::big_and_tiny_cube();
    auto result = meshseal::repair(mesh);
    // Tiny cube volume << 1% of big cube -> dropped
    CHECK(result.mesh.faces.size() < mesh.faces.size());
    CHECK(result.component_count == 1);
    bool dropped_note = false;
    for (const auto& n : result.notes)
        if (n.find("shells") != std::string::npos && n.find("dropped") != std::string::npos)
            dropped_note = true;
    CHECK(dropped_note);
}

TEST_CASE("repair: overlapping meshes resolved", "[e2e][intersections]") {
    auto mesh = meshseal::test::two_overlapping_tetrahedra();
    auto result = meshseal::repair(mesh);
    // Result should not be partial failure
    // (manifold should handle two separate tets even if they overlap)
    CHECK_FALSE(result.partial_failure);
    // Stage should have been applied
    bool has_intersect = false;
    for (const auto& s : result.stages_applied)
        if (s == "intersections") has_intersect = true;
    CHECK(has_intersect);
    // Output must be a valid closed solid.
    CHECK(boundary_edge_count(result.mesh) == 0u);
    CHECK(is_consistently_wound(result.mesh));
    CHECK(signed_volume(result.mesh) > 0.0);
}

TEST_CASE("repair: back-to-back triangles removed", "[e2e][thin]") {
    auto mesh = meshseal::test::back_to_back_triangles();
    // The two faces sit at z=0 and z=0.0005 with no shared vertex indices.
    // weld_tolerance=0.001 merges those near-coincident positions so both faces
    // share all three vertices and form one NO_BOUNDARY component (all edges valence 2).
    // degenerate=false prevents remove_duplicate_faces from eliminating the
    // antiparallel pair (canonical sort {0,1,2} matches for both faces) before
    // thin_features gets to inspect them.
    // orient=false preserves the original antiparallel windings.
    meshseal::RepairOptions opts;
    opts.weld_tolerance   = 0.001;  // merge z=0 and z=0.0005 vertices
    opts.degenerate       = false;  // prevent antiparallel pair being deduplicated
    opts.thin_features    = true;
    opts.orient           = false;  // preserve original antiparallel windings
    opts.holes            = false;  // don't fill holes (would seal the flat surface)
    opts.shells           = false;  // don't drop shells by volume
    opts.intersections    = false;  // don't run manifold self-union
    opts.soup_reconstruct = false;
    auto result = meshseal::repair(mesh, opts);
    // After repair, both faces should be removed
    CHECK(result.mesh.faces.size() == 0);
    // Stage applied
    bool thin_applied = false;
    for (const auto& s : result.stages_applied)
        if (s == "thin_features") thin_applied = true;
    CHECK(thin_applied);
}

TEST_CASE("repair: normal cube not affected by thin removal", "[e2e][thin]") {
    meshseal::RepairOptions opts;
    opts.thin_features = true;
    auto mesh = meshseal::test::unit_cube();
    auto result = meshseal::repair(mesh, opts);
    // Cube has no back-to-back faces within 1e-3 distance
    CHECK(result.mesh.faces.size() == mesh.faces.size());
    CHECK(result.watertight);
}

TEST_CASE("repair: cube soup reconstructed to manifold", "[e2e][soup]") {
    auto mesh = meshseal::test::cube_triangle_soup();
    // In the new component-aware pipeline the weld stage (default: on) merges the
    // 36 unshared vertex positions into 8 unique ones, turning the triangle soup
    // into a standard shared-vertex cube.  classify_components then sees one
    // 12-face NO_BOUNDARY component and routes it directly through the CLOSED
    // queue — soup_reconstruct is not needed (was_needed=false).
    // The key invariant: the output must be a valid closed solid.
    auto result = repair(mesh);  // all defaults
    CHECK_FALSE(result.partial_failure);
    REQUIRE(!result.mesh.faces.empty());
    CHECK(result.watertight);
    CHECK(result.is_volume);
    CHECK(boundary_edge_count(result.mesh) == 0u);
    CHECK(is_consistently_wound(result.mesh));
    CHECK(signed_volume(result.mesh) > 0.5);
}
