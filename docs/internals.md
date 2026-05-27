# meshseal — internals reference

Reference for the utility data structures and helpers under
`meshseal/src/internal/` and the I/O modules. Audience: someone reading
stage code wanting to understand "what does this helper actually do".
Pairs with `algorithms.md` (per-stage explanations).

## Mesh data model

The public type (in `include/meshseal/meshseal.h`):

```cpp
struct Mesh {
    std::vector<std::array<double, 3>> vertices;
    std::vector<std::array<uint32_t, 3>> faces;
};
```

Conventions:

- Vertices in **double precision**. STL/3MF inputs are float32 but the
  pipeline runs in double to preserve robust geometric predicates.
- Faces are **CCW outward** (positive signed volume when watertight).
- Face vertex indices are `uint32_t`. `uint32_t(-1)` is the conventional
  "no-such-face" / "INVALID" sentinel.
- No per-vertex / per-face attributes (color, UV) — meshseal is
  topology-only.

A float32 variant (`Meshf`) is provided for STL I/O round-tripping.

## I/O

### Binary STL (`meshseal/src/stl_io.cpp`)

The classic 84-byte header layout:

```
[0..80)   ignored header (80 bytes)
[80..84)  triangle count, uint32_t little-endian
[84..)    50 bytes per triangle: 3 floats normal + 9 floats verts + 2 bytes attr
```

The reader:

- Caps triangle count at **100M** to guard against malicious / corrupt
  count fields claiming `0xFFFFFFFF`.
- Validates `tri_count * 50 + 84 == file_size`; mismatch routes to ASCII
  detection (`looks_like_ascii_stl`).
- Reads each triangle's 9 float vert coords (the 3 normal floats are
  ignored — recomputed downstream from CCW winding).
- Per-face attribute bytes are read into `Mesh::attr` (used for color
  detection on multi-color exports; most readers ignore them).

The writer emits binary STL exclusively (ASCII STL is read-only).

### ASCII STL (`stl_io.cpp::read_stl_ascii`)

Some authoring tools export ASCII STL — `solid name … endsolid name`
with `facet normal x y z … outer loop vertex x y z … endloop endfacet`
records. The reader is lenient: it tokenises and grabs every `vertex
x y z` triple.

- Detection via `looks_like_ascii_stl`: file starts with `"solid"` AND
  size doesn't match the binary `84 + 50·count` layout (or `facet`
  keyword present in first 512 bytes).
- Hard cap of 300M vertex records (= 100M triangles) to match the
  binary path's DoS guard.
- Rejects NaN/Inf vertex coords with a clear error.

### 3MF (`meshseal/src/3mf_io.cpp`)

3MF is a ZIP archive containing:

- `[Content_Types].xml`, `_rels/.rels` — OPC packaging boilerplate.
- `3D/3dmodel.model` — the actual mesh data, XML.
- `Metadata/Slic3r_PE_model.config` — PrusaSlicer-specific per-volume
  metadata (if present): `<volume firstid="N" lastid="M" name="..."
  extruder="N">` blocks.

ZIP reading uses **miniz** (vendored, pinned to GIT_TAG `3.0.2` via
FetchContent — see CMakeLists.txt for the version-string idiosyncrasy
note).

**Zip-bomb guard.** Before `mz_zip_reader_extract_to_heap`, the reader
calls `mz_zip_reader_file_stat` and rejects entries whose declared
uncompressed size exceeds:

- 256 MB for the model XML (`3D/3dmodel.model`)
- 1 MB for the slicer config

This prevents a 1 KB malicious 3MF from advertising a 10 GB
uncompressed entry and forcing `malloc(10 GB)`.

**Locale-safe parsing.** All numeric attributes (`x`, `y`, `z`, `v1`,
`v2`, `v3`, `firstid`, `lastid`, `extruder`) use `std::from_chars`. The
former `std::stod` / `std::stoul` / `std::stoi` are locale-dependent;
on Czech / German Windows the decimal separator is `,` and `std::stod`
silently truncates `"1.5"` to `1`. The 3MF spec mandates `.`; `from_chars`
is the right tool.

**Parser caps.** 100M vertices and 100M triangles per `<mesh>` block.
NaN/Inf coords rejected.

**XML approach.** Hand-written substring scan (no proper XML parser).
Attributes extracted via `get_attr(tag, name)` which looks for ` name="`
(with leading space to avoid `fv1` matching `v1`). Limitation: legal
XML using tab/newline as attribute separator isn't accepted; not seen
in real-world 3MFs.

**Multi-volume 3MFs.** PrusaSlicer / Bambu Studio pack everything into
ONE `<object>` and put the partition into `Metadata/Slic3r_PE_model.config`.
The `read_3mf_volumes` API reads this metadata and returns each volume
as a separate sub-mesh (vertices re-indexed). The CLI uses this for
per-volume repair → per-`<object>` 3MF write (preserving multi-volume
structure round-trip).

The writer emits multi-`<object>` 3MF when given a vector of volumes,
with Slic3r config metadata for PrusaSlicer round-trip compatibility.

## SpatialHash (`internal/spatial_hash.h`)

Maps quantized (ix, iy, iz) cell coordinates to a list of indices.
Used by `weld`, `strip_doubled_membrane`, and other proximity-based
stages for O(1)-average lookup.

```cpp
class SpatialHash {
    explicit SpatialHash(double cell_size);
    void insert(const std::array<double,3>& p, uint32_t idx);
    std::vector<uint32_t> query_neighbors(const std::array<double,3>& p);
};
```

`query_neighbors` returns all indices in the 27-cell (3×3×3)
neighbourhood — guarantees coverage of any point within `cell_size` of
the query.

**Cell key**: a `CellKey { int64_t x, y, z; }` struct + custom hash
combiner (since 2026-05-27). The previous 21-bit-per-axis bit-packed
encoding wrapped every 2²¹ cells/axis, causing spurious bucket
collisions on extreme-aspect-ratio inputs (`lego/untitled2.stl`: 34 km
along X = 5 M cells / axis). Callers all metric-filter, so the wrap
was correctness-safe, but bloated buckets degraded performance.

## VoxelGrid (`internal/voxel_grid.{h,cpp}`)

Coarse 3D occupancy grid used as an inside/outside oracle by the
volumetric reconstructors (`stages/voxel_levelset`, `stages/alpha_wrap`,
`stages/soup_reconstruct`). Three operations:

1. **Rasterize** — mark cells overlapping each input triangle's AABB as
   `Surface`. Conservative (AABB, not exact Akenine-Möller SAT) but
   never misses a triangle.
2. **Flood-fill exterior** — 6-connected flood from every grid-boundary
   cell, marking reached cells as `Outside`. Remaining non-`Surface`
   cells are `Inside` (interior pockets).
3. **Sample** — return the label at a world point.

**Allocation safety**:

- Per-axis cap of 4096 cells.
- Per-axis sanity check `hi > lo`; degenerate axes get unit extent.
- Product cap of `2³⁰` = 1 GiB cells with halving fallback (halves the
  largest dimension and doubles cell_size until under cap).

The grid does NOT generate output geometry — output triangles come
from the original input, classified by sampling this grid (typically
at small offsets above and below their centroids).

## Generalized winding number (`internal/winding_number.{h,cpp}`)

Jacobson, Kavan, Sorkine-Hornung (SIGGRAPH 2013) robust formulation.

```cpp
double tri_solid_angle(a, b, c, p);  // steradians, signed
double generalized_winding_number(mesh, p);  // 0..1 for closed
```

The robust atan2 form (Van Oosterom & Strackee 1983):

```
Ω = 2·atan2(det[a-p, b-p, c-p],
            |a-p|·|b-p|·|c-p|
            + (a-p)·(b-p)·|c-p|
            + (b-p)·(c-p)·|a-p|
            + (c-p)·(a-p)·|b-p|)
```

For a closed orientable surface, ω(p) is `+1` inside, `0` outside. For
open / soup / multi-shell input, ω is fractional; the 0.5 level set is
a reasonable "best-guess inside" boundary.

**Cost**: O(F) per query. The orient_wn stage uses this per face → O(F²)
overall — acceptable up to ~50k F. For larger inputs a BVH-accelerated
hierarchical evaluation (Barill et al. 2018) would be needed; deferred.

**Defensive bounds check** (added 2026-05-26): silently skips faces
whose vertex indices exceed `vertices.size()` so a tainted face doesn't
silently UB across every query.

## arrangement (`internal/arrangement.{h,cpp}`)

Triangle-triangle intersection + co-linear cut consolidation +
triangle-by-cut-segment splitting. Used by `soup_reconstruct` Step 2
to discover natural meeting lines between overlapping triangles.

**Algorithm**: Möller-Trumbore for the intersection segment, then
`split_triangle_by_cuts` does the topological split.

**Scope of v1**:

- Non-coplanar pairs only. Coplanar pairs return no intersection. (The
  deferred coplanar_overlap plan would address this with constrained
  Delaunay re-triangulation.)
- The split assumes cut endpoints lie ON the triangle's edges (with
  tolerance). Cuts ending in the triangle's interior are dropped —
  would require a full PSLG / CDT.

**Tolerance caveats**: `kPlaneEps = kIntervalEps = 1e-12` are absolute,
compared against cross-products whose magnitude scales with `bbox²`. For
sub-mm CAD the test can pass spuriously; flagged in the review checklist
as a future scale-relative fix.

## self-intersection counter (`internal/si_count.{h,cpp}`)

Counts distinct triangle pairs that intersect, via Möller T-T on
spatial-hash-bucketed candidate pairs.

**Skip rule**: pairs sharing one or more vertices are NOT counted
(fan-around-vertex and edge-adjacent pairs produce numerical-artifact
"intersection segments" at the shared vertex/edge — not genuine SIs).
Matches the Python `si_count.py` convention.

**Caps** (for sweep timing):

- F ≤ 40 000: count is computed. Otherwise `-1` (= "not measured").
- Hit count capped at 10 000 pairs (once you're past that, the user
  already knows there's a problem).

## ManifoldWeld (`internal/manifold_weld.h`)

Welds a `manifold::MeshGL` output by **exact float32 bit-pattern key**.
Used by FWN-levelset, voxel_levelset, alpha_wrap, intersections.

Why exact float32: a manifold CSG output has float32 vertex coordinates
that may differ by 1 ULP between two faces that the CSG arithmetic
considered identical. A metric weld with tolerance `bbox·1e-7` would
either over-merge (too loose) or miss them (too tight). Exact bit
pattern matches the actual STL-format precision the downstream slicer
will see.

**`-0.0` normalisation** (added 2026-05-26): `+0.0f` and `-0.0f` have
different bit patterns. Before hashing, the weld now coerces signed-zero
to `+0`, so two sources of the same vertex with different zero-signs
dedup correctly.

## halfedge (`internal/halfedge.{h,cpp}`)

Half-edge mesh structure. Currently **dormant infrastructure** — no
stage in the current pipeline consumes it. Kept for future use (the
modular_pipeline_design might leverage it).

The build pre-condition is "no duplicate directed edges". To detect
violations, the build counts `directed_edge_collisions`; any future
caller can inspect this and reject malformed inputs rather than
silently keeping only the first-writer's half-edge.

## component_classifier (`internal/component_classifier.{h,cpp}`)

Edge-BFS components + open-ratio + planarity scoring + SOUP / OPEN /
NO_BOUNDARY classification used in Phase 5R.

**Open ratio** = (boundary edges) / (boundary + interior + NM edges).
Pure soup: ~1.0. Closed solid: 0.0.

**Planarity score**: relative thickness in the bbox's shortest axis —
near-zero for planar shells, near-one for solids.

The classifier is the gatekeeper for the heavyweight reconstructors;
its thresholds (`soup_open_ratio_threshold = 0.95`,
`soup_planarity_threshold = 0.01`) are deliberately conservative to
keep ordinary open-mesh repair on its own path.

## diagnostics (`internal/diagnostics.{h,cpp}`)

`compute_diagnostics(mesh)` returns a `MeshDiagnostics` with
`open_boundary_edges`, `non_manifold_edges`, `component_count`,
`signed_volume`, `is_orientable`, etc. O(F log F).

Called dozens of times per `repair()` invocation — every stage that
takes a defect-count or volume guard re-runs diagnostics on the pre-
and post-state.

A `cached_diag(mesh)` wrapper in `repair.cpp` memoises by mesh-pointer
+ face-count fingerprint to avoid repeated full recomputation on the
same mesh.

## vec3 (`internal/vec3.h`)

Plain 3-component double vector + dot, cross, length, sub, add, scale.
Header-only. Used everywhere; no `std::array` boilerplate at the use
site.

## timer (`internal/timer.h`)

`ScopedTimer` measures wall time of a stage and pushes the result into
`RepairResult::stage_times_ms`. RAII — the timer pushes on destruction.
Used by every stage call in `repair.cpp`.

The `MESHSEAL_PROFILE=1` env var enables additional `[prof]` per-stage
lines to stderr; zero overhead when unset.

## stage convention

Every stage under `meshseal/src/stages/` follows the same pattern:

```cpp
// stage.h
#pragma once
#include "../../include/meshseal/meshseal.h"
namespace meshseal::stages {

struct StageResult {
    Mesh     mesh;
    bool     applied      = false;
    uint32_t faces_changed = 0;
    // … per-stage metrics …
};

StageResult stage_name(const Mesh& mesh, /* optional params */);

} // namespace meshseal::stages
```

The orchestrator (`repair.cpp`) decides whether to adopt the result
based on the per-stage guard (typically: total-defect strictly
decreased, signed-volume preserved within some fraction).

**Adding a new stage**: see `CONTRIBUTING.md`-style notes in `CLAUDE.md`
under "Build quirks" — both `meshseal/CMakeLists.txt` AND the manually-
maintained `library/build/meshseal/meshseal.vcxproj` (gitignored) need
the new `.cpp` listed, because the local build tree's VS solution was
generated with a CMake generator the installed VS doesn't recognise.

## Numerical conventions

| concern | convention |
|---|---|
| **Internal precision** | double throughout |
| **I/O precision** | float32 (STL writes float32; 3MF same) |
| **Float32 round-trip** | post-pipeline snap+weld guarantees output survives float32 conversion |
| **Default tolerances** | scale-aware via `bbox_diag` or `mean_edge` multipliers |
| **Absolute floors** | `1e-9`, `1e-12`, `1e-30` for tolerance lower bounds |
| **Antiparallel test** | `n_a · n_b < −0.85` (default) or `< −0.95` (strict, for fold collapse) |
| **Coplanar test** | `|n_a · n_b| > 0.99` |
| **Vol-0 threshold** | `|signed_volume| < 1e-12` (replaced `bbox³ · 1e-6` which mis-classified thin parts) |

## Build system

`library/CMakeLists.txt` fetches three external dependencies via
`FetchContent`:

- **Catch2** v3.7.1 — unit testing.
- **manifold** v3.4.1 — CSG (used by `intersections`, `alpha_wrap`,
  `fwn_levelset`, `voxel_levelset`).
- **miniz** 3.0.2 — ZIP reading for 3MF.

**Build quirk** (CLAUDE.md "Build quirks"): the CMake-generated VS
solution uses generator `Visual Studio 18 2026` which the installed VS
Installer doesn't recognise. Do NOT regenerate. Build via MSBuild
directly:

```bash
cd library/build
MSBuild.exe meshseal/meshseal.vcxproj //p:Configuration=Release
MSBuild.exe cli/meshseal_cli.vcxproj //p:Configuration=Release
MSBuild.exe tests/tests_meshseal.vcxproj //p:Configuration=Release
./tests/Release/tests_meshseal.exe --reporter=compact
```

Forward-slash flags (`//p:...`) are required from Git Bash.

## Profiling

`MESHSEAL_PROFILE=1` env var → `[prof]` per-stage wall-clock + `[cnm]`
per-100-iteration progress logs inside `collapse_nm` to stderr.

`MESHSEAL_DUMP_DIR=<dir>` env var → per-stage STL dumps
(`dump_NN_<label>.stl`) inside `repair()`.

Both zero overhead when unset.

## Diagnostic infrastructure (external)

`C:/tmp/meshseal_diag/` (NOT in the repo) hosts Python tooling:

- `analyze.py` — STL + 3MF parser, per-mesh diagnostics
- `sweep.py` — runs the CLI on all corpus fixtures, tabulates CLEAN
- `surface_distance.py` / `surface_dist_sweep.py` — Hausdorff +
  mean-symmetric distance vs oracle
- `si_count.py` — Möller T-T self-intersection counter

These are dev apparatus, not shipped artifacts.
