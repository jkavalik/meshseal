#pragma once

#include <array>
#include <cstdint>
#include <map>
#include <optional>
#include <string>
#include <vector>

namespace meshseal {

struct Mesh {
    std::vector<std::array<double, 3>> vertices;
    std::vector<std::array<uint32_t, 3>> faces;
};

struct Meshf {
    std::vector<std::array<float, 3>> vertices;
    std::vector<std::array<uint32_t, 3>> faces;

    Mesh to_double() const;
    static Meshf from_double(const Mesh& m);
};

enum class ComponentStatus {
    repaired,
    unchanged_blocked,
    dropped_zero_volume,
    dropped_contained,
    failed_preserved,
    needs_reconstruction,
};

struct ComponentOutcome {
    int component_index;
    ComponentStatus status;
    std::string stage;   // stage that produced this outcome
    std::string reason;  // human-readable explanation
};

enum class EventSeverity {
    info,
    warning,
    error,
};

// Structured per-event record. Coexists with the free-form `notes` field
// for backward compatibility; new code that wants machine-readable
// diagnostics should prefer `events`. `code` is a stable short identifier
// (e.g. "weld.count", "holes.filled", "pile.preserved"), `message` is the
// same human-readable text that lands in `notes`.
struct RepairEvent {
    std::string    stage;
    EventSeverity  severity = EventSeverity::info;
    std::string    code;
    std::string    message;
};

struct RepairResult {
    Mesh mesh;
    bool watertight = false;
    bool is_volume = false;
    int component_count = 0;
    int self_intersections = 0;   // count of distinct triangle pairs that
                                  // self-intersect (Möller); -1 if not measured
    std::vector<std::string> stages_applied;
    std::map<std::string, double> stage_times_ms;
    std::vector<std::string> notes;            // free-form, backward compat
    std::vector<RepairEvent> events;           // structured equivalent of notes
    std::vector<ComponentOutcome> component_outcomes;
    bool partial_failure = false;
    // Confidence in the repair, in [0, 1]. Heuristic, not calibrated. Treat
    // as ordinal: 1.0 means "fully clean watertight solid", values < 0.5
    // mean the output is suspect (open boundaries, non-manifold edges, or
    // negligible volume). 0.0 means empty output.
    double confidence = 0.0;
};

struct RepairOptions {
    bool weld            = true;
    bool degenerate      = true;
    bool orient          = true;
    bool non_manifold    = true;
    bool holes           = true;
    bool shells          = true;
    bool intersections   = true;
    bool thin_features   = false;
    bool soup_reconstruct = true;

    // nullopt = auto-adaptive from mesh scale (bbox + mean edge length)
    std::optional<double> weld_tolerance;

    // Component classification thresholds (Phase 5R).
    // A component is routed to SOUP only if it is very heavily fragmented
    // (open_ratio close to 1.0) — either by the strict open_ratio==1 rule, or
    // by the planarity+ratio combination here. A loose 0.50 threshold caused
    // legitimate elongated open shells (e.g. an uncapped cylinder, with
    // open_ratio ~0.7 and small planarity) to be misclassified as SOUP and
    // sent to a CSG repair that cannot close them; 0.95 keeps the planarity
    // rule available for genuinely shattered planar inputs without
    // hijacking everyday open-mesh repair.
    double soup_planarity_threshold = 0.01;
    double soup_open_ratio_threshold = 0.95;

    // Internal: recursion depth for stages that re-enter repair() on a
    // locally-modified mesh (currently: nm_carve_refill, which carves the
    // NM-incident region and needs the full pipeline to settle the carved
    // mesh — Liepa + collapse_nm alone don't suffice on real-world inputs).
    // Bumped by the dispatcher; not for user code to set.
    int  recursion_depth      = 0;
    bool allow_carve_refill   = true;
};

RepairResult repair(const Mesh& mesh, const RepairOptions& opts = RepairOptions{});

} // namespace meshseal
