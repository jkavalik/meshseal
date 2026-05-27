# Changelog

All notable changes to meshseal are documented here. Format loosely
follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/);
versioning is [Semantic Versioning](https://semver.org/) — though the
**0.x** series makes no API-stability promises.

## 0.1.0 — 2026-05-27

First tagged release. Alpha quality; APIs not stable.

### Added

- **Repair pipeline** (~30 stages, orchestrated by `meshseal::repair`):
  input sanity validation, doubled-membrane strip, snap+weld,
  degenerate-face / duplicate removal, BFS winding-propagation orient,
  non-manifold vertex split + edge resolution, cross-shell boundary-
  loop bridging, component classification (SOUP / OPEN / NO_BOUNDARY),
  per-class handler routing, soup_reconstruct via Jacobson FWN or
  voxel-occupancy fallback, manifold CSG self-union, shell filtering
  by volume + containment, post-cleanup loop, T-junction splitting,
  NM-patch local remesh, progressive guarded edge-collapse,
  NM-local proximity weld, late thin-feature pass, junk-component
  drop, late targeted boundary-vertex merge, coplanar-fan drop,
  winding-number-guided orientation flip, NM carve + recursive
  repair() re-entry, alpha-wrap destruction fallback.
- **Mesh I/O**: binary + ASCII STL reader/writer; 3MF reader/writer
  with multi-`<object>` support and PrusaSlicer
  `Slic3r_PE_model.config` volume-partition round-trip.
- **Public API** (`<meshseal/meshseal.h>`):
  - `meshseal::Mesh` / `Meshf` data types
  - `meshseal::repair(mesh, opts)` entry point
  - `RepairOptions` with stage toggles, weld tolerance override,
    component-classification thresholds, recursion controls
  - `RepairResult` with diagnostics, stage timings, structured events,
    per-component outcomes, confidence score
  - `RepairOptions::on_progress` callback (per-stage events +
    cooperative cancellation)
- **CLI** (`meshseal_cli`): single-binary front-end, output format
  inferred from destination extension (`.stl` / `.3mf`), multi-volume
  3MF dispatch.

### Hardened against malicious input

- Phase 0 input sanity: NaN/Inf reject, bbox > 1e9 reject,
  face-index-out-of-range reject.
- Binary STL: 100M triangle cap.
- ASCII STL: 300M vertex-record cap.
- 3MF: 256 MB cap on uncompressed `3dmodel.model`, 1 MB cap on
  slicer config, 100M vertex + 100M triangle parser caps, NaN/Inf
  reject on parsed coordinates, locale-independent numeric parsing
  via `std::from_chars` / `fast_float`.
- All recursive `repair()` re-entries capped at depth 2.
- `std::manifold::Manifold` construction wrapped in try/catch.

### Cross-platform

- Builds on Windows (MSVC), Linux (GCC 10+, libstdc++), macOS arm64
  (Apple Clang, libc++). CI matrix exercises all three plus an
  older-GCC canary (GCC 10) verifying the `fast_float` polyfill path.
- `fast_float` polyfills `std::from_chars<double>` on toolchains where
  it's missing or `= delete`d (Apple Clang's libc++, libstdc++ < 11).
- STL binary I/O is explicit little-endian (works on big-endian hosts
  too — read/write helpers use byte shifts, not host-endian memcpy).
- Static-linked dependencies — `meshseal_cli` is a self-contained
  binary with no runtime `.dll` / `.so` requirements beyond
  libstdc++/libc++/libc.

### Quality

- **140 of 148 corpus fixtures CLEAN**; the remaining 8 are
  inherently unrepairable (vol-0 inputs, non-orientable surfaces,
  degenerate test probes).
- **634 + 36 unit-test assertions** across 116 + 4 test cases.
- **vs other open-source repair tools** (full 148-fixture corpus,
  default settings): meshseal 140 CLEAN, MeshFix 124, CGAL 111,
  pymeshlab 85, trimesh 77. See README for details.

### Known limitations

- **Linux determinism gap.** One corpus fixture (`black_vase`) lands
  CLEAN on MSVC but `nm=1` on libstdc++. The downstream
  `nm_carve_refill` rescue is calibrated to MSVC's
  `std::unordered_map` bucket iteration order; libstdc++ produces
  different greedy pair choices upstream. Documented; deferred.
- **Antiparallel manifold-edge pairs**: a small class of internal
  doubled-membrane defects with mismatched triangulations across the
  two sheets cannot be resolved by face-level surgery alone — they
  need a 2D CSG / constrained-Delaunay kernel (the deferred
  `coplanar_dedup` plan).
- **Intel-Mac coverage is not in CI** — GitHub's `macos-13` runner is
  effectively deprecated. macOS arm64 + Linux GCC 10 + Windows MSVC
  catch most issues a similar matrix would.
