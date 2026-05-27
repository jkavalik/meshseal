// Tests for RepairOptions::on_progress — the progress + cancellation
// callback wired through every prof_lap in repair().

#include "helpers.h"
#include <catch2/catch_test_macros.hpp>
#include <meshseal/meshseal.h>

#include <atomic>
#include <vector>

using namespace meshseal;

TEST_CASE("Progress callback: fires at every stage boundary", "[progress]") {
    auto mesh = test::unit_cube();
    RepairOptions opts;
    std::vector<ProgressEvent> events;
    opts.on_progress = [&](const ProgressEvent& ev) {
        events.push_back(ev);
        return true;  // continue
    };

    auto result = repair(mesh, opts);

    // At minimum, the final "done" event must have fired (a clean unit
    // cube triggers few prof_lap stages but always ends with "done").
    REQUIRE(!events.empty());
    CHECK(events.back().stage_name == "done");

    // elapsed_ms must be monotonic (or at least non-decreasing).
    for (size_t i = 1; i < events.size(); ++i) {
        CHECK(events[i].elapsed_ms >= events[i - 1].elapsed_ms);
    }

    // face_count of the final event should match result.mesh.
    CHECK(events.back().face_count == result.mesh.faces.size());

    // recursion_depth on the outermost call must be 0.
    CHECK(events.front().recursion_depth == 0);

    // Without cancellation, the repair must NOT be marked partial_failure
    // (a clean unit cube is trivially watertight).
    CHECK_FALSE(result.partial_failure);
    CHECK(result.watertight);
}

TEST_CASE("Progress callback: returning false cancels the pipeline", "[progress]") {
    auto mesh = test::unit_cube();
    RepairOptions opts;
    int seen = 0;
    opts.on_progress = [&](const ProgressEvent&) {
        ++seen;
        // Cancel on the first call.
        return false;
    };

    auto result = repair(mesh, opts);

    // Callback was called at least once before the cancel propagated.
    CHECK(seen >= 1);

    // Cancellation is reported as partial_failure with a recognizable note.
    CHECK(result.partial_failure);
    bool found_cancel_note = false;
    for (const auto& n : result.notes) {
        if (n.find("canceled") != std::string::npos) {
            found_cancel_note = true;
            break;
        }
    }
    CHECK(found_cancel_note);
}

TEST_CASE("Progress callback: default (empty) is a no-op", "[progress]") {
    auto mesh = test::unit_cube();
    RepairOptions opts;
    // opts.on_progress left as default (empty std::function)
    auto result = repair(mesh, opts);
    // Same expected outcome as the no-callback case.
    CHECK_FALSE(result.partial_failure);
    CHECK(result.watertight);
}

TEST_CASE("Progress callback: cancellation works on multi-volume input", "[progress]") {
    // Two separate cubes — exercises the spatial_split per-cluster
    // recursion path. Cancellation must propagate through inner repair()
    // frames cleanly.
    auto mesh = test::two_separate_cubes();
    RepairOptions opts;
    std::atomic<int> events_before_cancel{0};
    constexpr int kCancelAfter = 3;
    opts.on_progress = [&](const ProgressEvent&) {
        if (events_before_cancel.fetch_add(1) >= kCancelAfter) return false;
        return true;
    };

    auto result = repair(mesh, opts);

    CHECK(events_before_cancel >= kCancelAfter);
    CHECK(result.partial_failure);
}
