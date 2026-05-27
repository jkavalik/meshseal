# meshseal — algorithms reference

Stage-by-stage explanation of the repair pipeline. Audience: someone
reading the code wanting to understand "what is this stage doing and
why is it shaped this way". Pairs with `internals.md` (data structures)
and the project root `CLAUDE.md` (decision log + history).

The pipeline orchestrator is `meshseal/src/repair.cpp`. Each stage is a
self-contained translation unit under `meshseal/src/stages/` exposing a
single entry function that takes a `Mesh` (or `const Mesh&`) and returns
either a new `Mesh` or a small result struct.

The order matters. The pipeline is not a tree of independent passes —
each stage assumes the previous one left the mesh in a particular state.
The phase numbers below match the labels in `repair.cpp`.

## Mesh data model

```cpp
struct Mesh {
    std::vector<std::array<double, 3>> vertices;  // double-precision
    std::vector<std::array<uint32_t, 3>> faces;   // CCW outward winding
};
```

Faces reference vertices by index. Winding convention is CCW seen from
outside the solid (positive signed volume for a watertight mesh).

Throughout the pipeline:

- **bnd** = boundary edges, count of edges with exactly one incident face.
- **nm** = non-manifold edges, count of edges with ≥ 3 incident faces.
- **antipar** = manifold edges (exactly 2 incident faces) whose two
  normals are antiparallel (`n_a · n_b < -0.85`); a silent quality
  metric, not part of the strict CLEAN check.
- **CLEAN** = `bnd == 0 && nm == 0 && |vol| > 1e-12`.

---

## Phase 0 — input sanity

`repair.cpp:85–145`. Three rejects before any expensive work:

1. **NaN/Inf vertex coords** — corrupt files (real D-batch example:
   `apple/usb_end.stl`, 1e34-magnitude floats) used to spin the pipeline
   for minutes.
2. **`bbox_diag > 1e9`** — ~1000 km. Any larger is unit-scale corruption.
3. **Face index out-of-range** — library-API safety: a caller building
   a `Mesh` by hand could put indices past `vertices.size()`. Downstream
   `mesh.vertices[face[k]]` would be undefined behaviour.

On reject: `result.partial_failure = true`, empty mesh, descriptive note.

## Phase 0.5 — strip_doubled_membrane

`stages/strip_doubled_membrane.cpp`. Runs **before weld** because the
defect it targets relies on the two doubled face sheets having
distinct (but coincident-after-weld) vertex sets.

**What it removes.** Multi-color / multi-shell exports (PrusaSlicer,
Bambu Studio, MeshMixer) frequently contain an **internal doubled
membrane** — where two originally-separate solid bodies were left
touching during export, the contact interface survives as two
back-to-back face sheets with antiparallel normals. Without this stage,
weld merges the coincident vertex positions but the existing topology
pipeline only ever strips **one** of the two sheets, leaving the
survivor as an orphaned internal wall.

**Algorithm.** For each face F with centroid C and unit normal N, find
any face F' such that `||C_F − C_F'|| < tol·bbox_diag` AND
`N_F · N_F' < −0.85`. Greedy pairing: mark both as membrane if neither
is already marked. Spatial hash (cell = 2·tol) gives O(F) average query.

**Three gates, ALL required to adopt:**

- `total_faces ≥ 100` — tiny synthetic test meshes (a 2-face back-to-back
  pair, or the 24-face big-and-tiny-cube test fixture) would otherwise
  trigger the strip on legitimate small objects.
- `detected_fraction ≥ 0.10` — small localised patches (trumpet's 3 %
  annular seam, black_vase's 0 %) fall below the gate and continue
  through the existing `collapse_nm` / `nm_local_repair` pipeline. A
  wholesale strip on a 3 % patch creates boundary loops that Liepa
  fills in ways that introduce more antipar than they remove.
- `detected_fraction ≤ 0.90` — fully-doubled scan reconstructions
  (`room_surface_mesh.stl`, 100 % back-to-back from the scanner picking
  up both sides) skip the strip and route through `soup_reconstruct`
  instead; stripping 100 % leaves an empty mesh.

**Concrete win.** Bee_v3.stl: F=78920 nm=44 antipar=1074 in 21 s →
F=22990 nm=0 antipar=49 in 1.8 s (12× faster, 95 % antipar reduction).

## Phase 1 — weld

`stages/weld.cpp`. Merge vertices that share a position to within a
tolerance. Tolerance is **scale-aware** by default: the auto-tolerance
is `max(1e-9, bbox·1e-7, mean_edge·1e-4)`.

**Why the mean-edge term.** A pure `bbox·1e-7` tolerance is dwarfed by
geometric features on dense meshes (the mean edge length encodes the
feature size). `mean_edge·1e-4` ensures within-feature precision.

**Why not a tighter weld.** Real-world doubled-vertex bugs (e.g.
multi-color exports) sit at distances ~1e-7 of bbox; a tighter weld
misses them. A looser weld over-merges nearby legitimate vertices.

The weld uses a spatial hash (`internal/spatial_hash.h`) with cell size
= 2 · tolerance to guarantee 27-cell neighbour coverage.

## Phase 2 — degenerate

`stages/degenerate.cpp`. Removes two classes of degenerate face:

1. **Zero-area** — area below `bbox² · 1e-22` absolute threshold (the
   header default `1e-14` is an absolute fallback). Threshold of 0.0
   triggers scale-aware behaviour.
2. **Duplicate faces** — exact vertex-set match (irrespective of
   winding) → keeps one canonical copy. Catches back-to-back triangle
   pairs whose vertex triples are identical after weld.

Why "after weld": before weld, two coincident-position faces have
distinct vertex indices and aren't detected; the canonical-sort dedup
needs shared vertex indices.

## Phase 3 — orient

`stages/orient.cpp`. Propagates a consistent CCW winding across each
connected component.

**Algorithm.** BFS over face adjacency:

1. Build directed-edge map `(u, v) → face_index`.
2. Pick an unvisited seed face.
3. BFS: for each neighbour `nb` of the current face, check the shared
   edge's orientation. If `nb` has the SAME directed edge as current,
   their windings disagree → flip `nb`. If `nb` has the OPPOSITE
   directed edge, windings agree — leave alone.
4. After BFS settles each component, if the component's signed volume
   is negative, flip ALL faces in the component (per-component
   sign-flip; see `do_signed_volume_flip`).

The signed-volume flip is **destructive** for non-orientable inputs
(Möbius strip) because it interprets the un-orientable surface as
needing inversion. Phase 4c re-runs orient with `do_signed_volume_flip
= false` after `nm_edge` cleans up — propagates winding across newly
manifold-by-then edges without re-flipping.

## Phase 4 — non-manifold cleanup

Three substages in order:

### 4a — nm_vertex

`stages/nm_vertex.cpp`. Splits **non-manifold vertices** (where two or
more face fans meet at a point but share no edges). Each fan gets its
own copy of the vertex, separated by a tiny offset
(`1e-6 · fan_radius`) toward the fan centroid. Float32-stable: the
offset survives the float32 round-trip on slicer write.

### 4b — nm_edge

`stages/nm_edge.cpp`. Resolves **non-manifold edges** (3+ incident
faces). For each NM edge, separates faces into directed-edge groups:
group A traverses u→v, group B traverses v→u. Two cases:

- **Fin removal**: when one group is dominant (large coherent surface)
  and the other has stragglers (1-2 fin faces with normals
  incoherent with the dominant group), the stragglers are dropped.
- **Multi-surface preserve**: when BOTH groups have ≥ 2 faces AND
  intra-group normal coherence ≥ 0.85 AND fans ≤ 3 AND input ≥ 100
  faces, the NM edge is preserved as a genuine multi-surface meeting
  (e.g. three planes meeting at a T-edge — the edge legitimately has
  3+ incident faces).

The greedy resolution depends on `unordered_map` iteration order; this
non-determinism is documented and the downstream pipeline (especially
the `nm_carve_refill` rescue for black_vase) is calibrated to the
current bucket-order behaviour. Sorting was tried and reverted —
see `REVIEW_CHECKLIST.md` Phase 3.

### 4c — orient (BFS-only)

Re-run of `stages/orient.cpp` with `do_signed_volume_flip = false`.
Propagates winding across edges that became manifold after 4b.

### 4d — bridge_loops

`stages/bridge_loops.cpp`. Zips two coincident open boundary loops from
different components. Targets the "cap on a separate cavity shell"
pattern (captain_toad's mouth-cavity is a separate shell from the head;
the mouth opening is a boundary loop on each).

**Gates** (all required, all empirically necessary):

- Both loops from **different** connected components.
- Centroid-to-centroid distance `< 0.35 · average_radius`.
- Loop-radius ratio `≥ 0.5` (similar size).
- Vertex-count ratio `≥ 0.5` (similar discretization density).
- Loop normals' dot product `< −0.3` (facing each other).
- Möller line-of-sight unobstructed between centroids.

Plus a per-strip safety guard: roll back any strip that raises nm-edge
count or fails to lower boundary count.

The loose-gate Go-port version regressed bird_bath / human / sea_vase /
black_vase by zipping unrelated holes.

## Phase 5R — component classification

`internal/component_classifier.cpp`. Classifies each connected component
into one of three buckets:

- **SOUP** — heavily fragmented (open_ratio close to 1.0). Routed to
  `soup_reconstruct` / FWN-levelset / voxel-levelset rebuild.
- **OPEN** — has open boundaries but is recognisably one piece (low
  open_ratio, structure preserved). Routed to per-component nm_edge +
  orient + Liepa hole-fill.
- **NO_BOUNDARY** — already closed. Routed to manifold-CSG self-union.

The thresholds (`soup_open_ratio_threshold = 0.95`,
`soup_planarity_threshold = 0.01`) were chosen empirically — a loose
0.50 threshold caused legitimate elongated open shells (uncapped
cylinder, open_ratio ~0.7) to be misclassified as SOUP and sent to a
CSG repair that cannot close them.

## Phase 6R — soup_reconstruct (SOUP path)

`stages/soup_reconstruct.cpp`. Per-component reconstruction for
classified-SOUP components. Two reconstructors with a Stage-10 fallback:

- **FWN-first**: Jacobson 2013 generalized winding number computed at
  the centroid of each input face → fed to `manifold::Manifold::LevelSet`
  for marching-cubes extraction. Smoother surface than voxel; better
  for Phase 8R rescue.
- **Voxel-occupancy fallback**: when FWN's output is empty / degenerate,
  fall back to the orientation-agnostic voxel-occupancy SDF. This is
  the originals' prescribed algorithm for random-winding soup.

Both reconstructors share `manifold_to_welded_mesh()` (welds MeshGL
output by bit-exact float32 key — see `internals.md`).

The `fwn_used` flag gates downstream `intersections` off because
re-running manifold CSG on dense marching-cubes output would re-introduce
NM at float32 collisions.

## Phase 7R — OPEN handler

Per-component for classified-OPEN. Runs `nm_edge` (drop fins) + `orient`
+ `holes::fill_holes` (Liepa) on the component in isolation. Closed
components are pushed into the manifold-closed queue for Phase 9R.

## Phase 8R — closed_queue

For classified-NO_BOUNDARY plus rescued-OPEN components: `orient` +
`nm_edge`. Closed components feed into `manifold_closed_queue`;
unrescuable ones become `soup_fallback`; tiny leftovers go to `pile`.

**Catastrophic-collapse guard** — a real-input component of ≥ 1000
faces that emerges from `soup_reconstruct` with < 10 % of its input
face count is presumed destroyed and rejected (preserve original
instead). Distinguishes real input from FWN/soup output via the
`closed_queue_preserve` flag.

## Phase 9R — intersections

`stages/intersections.cpp`. Self-union of all manifold-ingestible
components via `manifold::Manifold::Decompose() + sequential union`.
Resolves CSG-style interpenetration.

**Gates**: `component_count > 1 || nm > 0`. Skipped if `fwn_used`
(see Phase 6R).

**AABB-disjoint pre-check** — when components are AABB-disjoint, the
union is provably a no-op. Skip with a note.

**Try/catch around manifold construction** — the library is documented
to throw on some pathological inputs; a throw returns `manifold_failed
= true` instead of propagating up through every outer `repair()` frame.

### Shells

`stages/shells.cpp`. Filters components by volume and containment:

- Drop components with `|vol| < vol_threshold·max_vol` (= 0.001·max_vol).
- Drop components fully contained in the main body (ray-cast +X, odd
  crossings); the threshold is `≤ max(200, nf/100) faces, |vol| < 1 %
  body`.

## Phase 10R — orient + degenerate (cleanup)

A pass of `orient_mesh` (without sign-flip) and `remove_degenerate` to
clean up post-CSG artifacts.

### 10R.5 — pile reintegration

Tiny components that didn't survive the OPEN/NO_BOUNDARY split are held
in a `pile`. If `soup_reconstruct` succeeds on the pile, the
reconstructed result is folded back in. Synthetic FWN-derived entries
are NOT re-added (they're MC garbage that would re-explode).

## Float32-compatibility pass

The output STL/3MF is written in float32. To ensure the float32 round-trip
doesn't re-introduce NM (two distinct double-precision verts that snap to
the same float32 value), the pipeline does a final `snap + weld +
optional manifold rebuild`.

## Post-cleanup loop (≤ 3 iterations)

Cycles `nm_edge → orient → fill_holes` while total defect strictly
decreases. Snapshot-and-revert if any iteration makes things worse.

Boundary threshold raised to 5000 in 2026-05 — previously 400, which
refused to even try on multi-color 3MF outputs landing just above the
old cap. The snapshot-revert guard inside the loop is the real safety
net.

## tjunction

`stages/tjunction.cpp`. Splits non-conforming collinear-vertex edges.

A T-junction is one face spanning edge (a, b) while the other side
splits it into sub-edges through vertices lying collinear on (a, b)'s
interior. It reads as a residual boundary or NM edge that hole-fill
cannot close (a zero-width collinear slit is a zero-area loop).

The stage scans problem edges (count ≠ 2) for collinear interior
vertices (perpendicular distance `< bbox·1e-6`, param strictly in
(0, 1)) and fan-splits each incident face through them: `(a, b, c) →
(a, v1, c), (v1, v2, c), …, (vk, b, c)` — winding-preserved, no
degenerate faces.

## nm_patch_remesh

`stages/nm_patch_remesh.cpp`. For residual non-manifold edges (≥ 3
incident faces) it deletes a small BFS-expanded patch of faces around
each NM edge and re-closes the opened boundary loops with Liepa
hole-fill → a clean single-layer 2-manifold disk.

**Critical prerequisite.** The residual NM on doubled-surface fixtures
is *index-space-invisible* — the doubled layer uses distinct vertex
indices whose positions are bit-identical only once truncated to
float32. So `compute_diagnostics` (index-based) reports nm=0 while
analyze (exact-float32 dedup) sees the NM. The block first snaps to
float32 and welds by exact float32 bit pattern (NOT a `bbox·1e-7`
metric weld — that misses them) so the NM becomes visible to the
index-based stages.

**Ring-retry**. Tries multiple BFS ring sizes per iteration
(`{1, 0, 2, 3, 4}`) and picks the smallest one that strictly improves
total defect.

**Two passes**. A second nm_patch_remesh pass runs after `collapse_nm`
— the collapse_nm-cleaned mesh sometimes admits patches the earlier
pass couldn't break into.

## collapse_nm

`stages/collapse_nm.cpp`. Progressive guarded edge-collapse for
residual NM flap geometry. Per-collapse guards: NM-edge count and
open-boundary count may only fall or hold (watertight preserved at
every step), no non-sliver face may flip normal.

**Two phases:**

- **Phase 1** erases NM edges. Two tiers: "improve" (NM strictly drops)
  and "neutral" (NM/bnd hold, region face count shrinks — decimate to
  escape a stall). Phase 1 has a stall-streak cap (`kStallCap = 50`) —
  break after 50 consecutive non-improving iterations.
- **Phase 2** erases residual zero-volume flaps once Phase 1 reaches
  `nm == 0`. A "flap" is a manifold edge with antiparallel adjacent
  faces (`n_a · n_b < −0.95`). Phase 2 collapses carry a near-exact-zero
  per-collapse volume cap (`|Δvol| ≤ 1e-9·total`) — a fold is erased
  only when it bounds provably zero volume.

**Performance knobs** (kytka1 perf cascade from 2026-05-23):

- Candidate gather is O(NM-region), not O(E). Iterates `v2f[x]` for
  `x ∈ nmv` only, not all edges.
- v2f-compaction trigger fires when stale-entry count > 2 · alive-entry
  count and > 1024 — rebuilds `v2f` from the alive face list.
- The `max_collapse_len = 0.006 · bbox_diag` cap. Must stay below wall
  thickness or it welds opposite walls of thin-wall models (a 0.05
  cap ballooned black_vase volume +52 %).

## nm_local_repair

`stages/nm_local_repair.cpp`. Two combined phases targeting NM edges only:

1. **NM-local proximity weld**: for each NM-edge endpoint, merge any
   1-ring neighbour within `bbox·1e-4`. Catches CSG-corner near-coincident
   duplicates the global weld (`bbox·1e-7`) is too tight to handle
   (`t10k_1582375` has chamfer corners 0.004 mm apart = 0.007 % of
   bbox — should be one vertex).
2. **Strict back-to-back dedup**: exact `frozenset(vertices)` match
   with opposite winding. Strictly identifies true duplicates;
   `remove_thin_features`' centroid+normal metric was rejected by the
   defect guard on `bumpy_white` for catching non-pair faces.

Combined defect+volume guard (≤ 0.5 % vol change).

## Late stages

Run unconditionally after collapse_nm, defect-guarded:

- **thin_features** (`stages/thin_features.cpp` via `remove_thin_features`) —
  back-to-back face pair removal by centroid + anti-parallel normal.
  Catches Liepa fill's accidental fold pairs.
- **junk-component drop** — edge-BFS into connected components, then
  drops: tiny near-zero-volume flaps (≤ 30 faces, `|vol| < 1e-5·max_vol`,
  back-to-back duplicate seam junk); small components fully CONTAINED
  in the main body (ray-cast, ≤ max(200, nf/100) faces, `|vol| < 1 %`
  body).
- **late_fill_holes** — targeted boundary-vertex merge for collapse_nm
  pinch residuals. Pairs of boundary vertices within `bbox · 2e-7`
  (~2 float32 ULPs) merge via face-index remapping. Restricted to the
  boundary-vertex set so the tight tolerance can't cause downstream
  over-merge.
- **coplanar_fan_drop** — drops same-sided coplanar 4-fan duplicate
  pairs at NM edges; limited-stage recursive cleanup.

## orient_wn

`stages/orient_wn.cpp`. GWN-based per-patch orientation flip. Gated on
`antipar_pre ≥ 50 AND nm == 0`; dormant on most corpus inputs.

**Algorithm.** Patch-BFS through consistent-winding manifold edges
(antipar edges are patch boundaries) + per-face relative-GWN majority
vote: for each face sample c+eps*N and c−eps*N; correctly-oriented
face has `GWN_out < GWN_in`. Flip patches whose majority of faces vote
flip.

**Four guards, ALL required to adopt:**

- `nm_pre == 0` — running on meshes with NM caused regressions on
  `kuzely.stl`.
- antipar count strict decrease.
- total signed-volume sign not flipped — Bee_v3 vol[3] tried to flip
  99 % of the main patch without it.
- NM edge count not increased.

## nm_carve_refill — recursive re-entry

`stages/nm_carve_refill.cpp` + the carve loop in `repair.cpp`. Handles
small stubborn residual NM by carving the NM-incident faces + 1-ring
halo, then re-entering `repair()` recursively on the carved mesh.

**Gates:**

- `opts.allow_carve_refill` true (recursive call sets false).
- `opts.recursion_depth < 2`.
- Position-visible NM ≤ 20 (gate against kytka1-class spirals — Liepa
  on > 50-NM input produces > 100-NM result that grinds for minutes).
- Recursive result strictly improves total defect AND stays within
  5 % volume (dual-scale: signed-volume OR bbox-volume, for thin parts
  where signed-volume is unreliable).

The recursion runs the FULL pipeline on the carved mesh — Liepa +
collapse_nm alone can't replicate it on real-world inputs.

## destruction fallback (alpha_wrap)

`stages/alpha_wrap.cpp`. Wholesale reconstruction for fragmented-solid
inputs (LDraw / LEGO soup — separate shells exported as one mesh).

**Algorithm.** Morphological closing with ball of radius α:

1. Voxelise (longest-edge-bisection rasterisation).
2. Distance transform #1: separable F&H EDT.
3. Flood-fill the reachable probe-centre set C through
   `{dist_to_surface > α}`.
4. Distance transform #2.
5. Continuous SDF `s = dist_to_C − α` (positive inside).
6. `manifold::LevelSet` extracts the iso-surface.
7. Exact-float32 weld.

Wired into the destruction fallback in `repair()`: if the pipeline
output is non-watertight AND fragmented (`component_count > 1`), run
alpha_wrap and adopt only if clean watertight + > 1.5× volume recovered.

The `component_count > 1` pre-gate is the key — a single-component
near-miss (`black_vase` nm=3, `t10k_*` nm=4 — full volume, good
geometry) is left untouched.

## do-no-harm guard

Final stage of `repair()`. Records at entry whether the input is
already a clean watertight solid (float32-welded, `bnd=0 && nm=0 &&
|vol| > 1e-12`); if so AND the pipeline ended with a non-watertight
result, restores the float32-welded input unchanged.

A repair tool must never hand back something worse than it got.
Real-world trigger: inputs that are edge-manifold and closed but have
non-manifold *vertices* (bowtie pinches) — the `nm_vertex` split, an
over-aggressive `weld`, or the `intersections` manifold self-union
introduce NM the later stages cannot fully clear.

## spatial_split (multi-object STL dispatcher)

`repair.cpp`. Before the main pipeline, if the input has AABB-disjoint
clusters of connected components (typical of a multi-object print
plate exported as one STL), recurse per-cluster with
`allow_spatial_split = false`.

**Triggers**: `clusters.size() >= 3` OR `clusters.size() >= 2 && F > 1M`.

**Why**: alpha_wrap voxel grid is sized per-cluster instead of
whole-plate. NM count is unchanged (AABB-disjoint shells have zero
geometric coupling). Single-cluster inputs (most of the corpus) fall
through at near-zero cost.

**post_split_carve_refill** — after merging cluster outputs, if the
merged result still has small residual NM (1 ≤ nm ≤ 20), run
`nm_carve_refill` once on the merged mesh and adopt the recursive
repair if it strictly improves. Concrete fix: snorlax.3mf (14-cluster
vol[2] where cluster 8 left nm=1) → CLEAN.

---

## Dormant infrastructure

Stages compiled in but not currently wired into `repair()`:

- **`sliver.cpp`** — extreme-needle triangle detector + shortest-edge
  collapse. Pillinger quality `q = 4√3·A/Σℓ²` below 0.02 triggers
  collapse of the shortest manifold edge. Kept for genuine needle-sliver
  inputs; the corpus's residual-NM cases are structural (doubled
  membranes), not sliver-shaped, so the stage is a no-op on every
  current fixture.
- **`si_split.cpp`** — standalone arrangement-based SI split (Möller
  intersection + co-linear cut consolidation + `split_triangle_by_cuts`).
  Doesn't include an inside/outside selection step, so used standalone
  it makes meshes worse (splitting creates new NM at the cut edges).
  Kept as future infrastructure if a selection step is built.
- **`fwn_levelset.cpp`** — Jacobson FWN driving `manifold::LevelSet`.
  Called by `soup_reconstruct` (Phase 6R) but also exported standalone
  for the destruction fallback's old code path. Now subsumed by
  `alpha_wrap`'s morphological-closing approach which is the active
  fallback.
- **`halfedge`** (internal) — half-edge data structure with no current
  callers. Has a `directed_edge_collisions` counter so a future caller
  can detect precondition violations.

## Algorithm references

- Jacobson 2013, "Robust Inside-Outside Segmentation using Generalized
  Winding Numbers" — Van Oosterom & Strackee solid-angle formula used
  in `internal/winding_number.cpp` and the FWN-levelset stage.
- Liepa 2003, "Filling Holes in Meshes" — the DP-based hole filler in
  `stages/holes.cpp`. Dihedral-weight term is implemented but disabled.
- Möller, "A Fast Triangle-Triangle Intersection Test" — used by
  `internal/arrangement.cpp` and `internal/si_count.cpp`.
- Felzenszwalb–Huttenlocher 2004, "Distance Transforms of Sampled
  Functions" — separable EDT in `stages/alpha_wrap.cpp`.

Project-level design docs (in repo root `docs/`) cover the
alpha-wrapping plan, the deferred coplanar_overlap (2D-CSG) plan, the
modular_pipeline_design, and historical algorithm investigations.
