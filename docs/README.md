# meshseal documentation

Reference docs for the meshseal mesh-repair library.

- **[`algorithms.md`](algorithms.md)** — per-stage explanation of the
  repair pipeline. What each stage does, why it's shaped the way it is,
  what it depends on from preceding stages.
- **[`internals.md`](internals.md)** — utility data structures,
  numerical conventions, I/O formats, build system.

For history, decision log, and the running design notes that the
project lives by, see `../../CLAUDE.md` in the repo root.

For design plans (alpha-wrapping, modular pipeline, deferred coplanar
overlap kernel, shell isolation), see `../../docs/` in the repo root.

## Quick orientation

If you've just opened the codebase:

1. Start with `algorithms.md` "Mesh data model" + "Phase 0 — input
   sanity" sections (5 minutes) to know what the pipeline operates on.
2. Skim the phase list in `algorithms.md` to get the pipeline shape.
3. Open `meshseal/src/repair.cpp` and read top-down, jumping to
   `algorithms.md` whenever a stage name is unfamiliar.
4. Open individual stages in `meshseal/src/stages/` as needed; each is
   self-contained.

For utility code (BVH, spatial hash, voxel grid, half-edge, winding
number, etc.) read `internals.md` first.
