# meshseal

A C++17 mesh repair library and CLI for 3D printing. Reads STL or 3MF,
outputs a watertight 2-manifold solid suitable for slicing.

> **Status: 0.1.0 — alpha.** APIs are not stable; the library is being
> developed against a real-world corpus of 148 defective fixtures plus
> ~1200 anomalous real-world STLs from a developer's print archive.
> Cross-platform CI on Windows, Linux, and macOS.

## What it does

Slicing algorithms (FDM, SLA) require a closed manifold mesh. Real-world
STL files routinely violate this: open boundaries, non-manifold edges,
duplicated face pairs, inverted normals, fragmented solids, and so on.
meshseal pipelines ~30 repair stages to clean these up while
**preserving the original triangulation wherever possible**.

The current corpus result: **140 of 148 fixtures CLEAN** (the remaining
8 are inherently unrepairable: vol-0 inputs, non-orientable surfaces,
degenerate test probes). On a real-world batch of ~424 anomalous STLs
from print archives: **~97% CLEAN** end-to-end.

See [How it compares](#how-it-compares) for benchmark numbers vs
other open-source repair tools.

## Quick start

### CLI

```sh
# Build (requires CMake 3.18+, C++17 compiler)
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --config Release

# Repair an STL
./build/cli/meshseal_cli   input.stl  output.stl
./build/cli/meshseal_cli   input.3mf  output.3mf       # 3MF in/out, multi-volume preserved
./build/cli/meshseal_cli   input.stl  output.3mf       # cross-format
```

Output format is chosen by the destination extension. For 3MF input
with multi-volume metadata (PrusaSlicer-style), each volume is repaired
independently and written back into its own `<object>` in the output 3MF.

### Library

Minimal example:

```cpp
#include <meshseal/meshseal.h>
#include <iostream>

int main() {
    meshseal::Mesh in;
    // ... populate in.vertices and in.faces from your source ...

    meshseal::RepairOptions opts;
    meshseal::RepairResult res = meshseal::repair(in, opts);

    if (res.watertight) {
        std::cout << "repaired: F=" << res.mesh.faces.size()
                  << "  confidence=" << res.confidence << "\n";
    } else {
        std::cout << "partial: ";
        for (const auto& n : res.notes) std::cout << "  " << n << "\n";
    }
}
```

CMake integration via FetchContent:

```cmake
include(FetchContent)
FetchContent_Declare(meshseal
    GIT_REPOSITORY https://github.com/jkavalik/meshseal.git
    GIT_TAG        v0.1.0)
FetchContent_MakeAvailable(meshseal)

target_link_libraries(my_app PRIVATE meshseal)
```

### Progress reporting and cancellation

For long-running repairs (large meshes, complex defect classes),
attach a callback:

```cpp
opts.on_progress = [](const meshseal::ProgressEvent& e) -> bool {
    std::cout << "[" << e.elapsed_ms << " ms] " << e.stage_name
              << "  F=" << e.face_count << "\n";
    return !user_clicked_cancel;  // return false to bail
};
```

Cancellation propagates through recursive sub-repairs and returns with
`partial_failure = true`.

## Supported platforms

CI is green on:

- **Windows** — MSVC 2022 (Visual Studio 17)
- **Linux** — GCC 13 (ubuntu-latest) and GCC 10 (ubuntu-22.04 canary)
- **macOS** — Apple Clang on arm64 (macos-latest)

Intel-Mac is not currently tested in CI — GitHub deprecated the
free `macos-13` runner. The library should work on Intel-Mac; user
reports needed to confirm.

Minimum requirements: C++17 compiler, CMake 3.18+. No external system
libraries beyond libstdc++ / libc++ / libc.

## How it compares

Run on the full 148-fixture corpus across five mesh-repair tools with
their default settings (Linux WSL, 60 s per-tool timeout):

| | meshseal | meshfix | cgal | pymeshlab | trimesh |
|---|---|---|---|---|---|
| **CLEAN** | **140 (94.6%)** | 124 (83.8%) | 111 (75.0%) | 85 (57.4%) | 77 (52.0%) |
| Total time | **47 s** | 78 s | 82 s | 14 s | 17 s |

- **MeshFix** ([Attene 2010](https://doi.org/10.1007/s00371-010-0416-3))
  is the closest competitor — same geometry-preserving design. meshseal
  leads by 16 fixtures and is also faster.
- **CGAL** is a general geometry library; we ran the default
  `repair_polygon_soup + orient + triangulate_holes` pipeline. A
  custom pipeline with aggressive non-manifold-edge removal would close
  some of the gap.
- **pymeshlab / trimesh** are lightweight Python wrappers — fast but
  much less aggressive on the harder defect classes.
- The standout category is **triangle_soup: 15/15 vs 2/15** for the
  other three serious tools — meshseal's reconstruct stage (FWN +
  alpha-wrap) handles fragmented input the others don't try.
- On **real-world** fixtures: meshseal 36/37, MeshFix 35/37,
  CGAL 32/37 — competitive at the top.

Caveats: all tools at default settings (skilled users could tune any
of them further); CLEAN check requires `bnd=0 && nm=0 && |vol|>1e-12`
— some lightweight tools "win" CLEAN on pathological inputs by
collapsing them to a unit tetrahedron, which meshseal refuses
(geometry-preservation principle). One `black_vase` regression on
Linux is the documented `unordered_map` iteration-order gap; meshseal
on Win/Mac scores 141/148.

## What's in v0.1.0

See [`CHANGELOG.md`](CHANGELOG.md) for the full release notes.

Highlights:

- ~30-stage repair pipeline with conservative geometry-preservation
  guards: snap+weld → degenerate-removal → orient → non-manifold
  cleanup → component classification → soup reconstruction → CSG
  union → shell filtering → late iterative cleanup → carve+refill
  recursion → alpha-wrap fallback.
- Multi-volume 3MF round-trip (reads & writes PrusaSlicer's
  `Slic3r_PE_model.config` volume metadata).
- Progress + cancellation callback API.
- ASCII + binary STL reader/writer, 3MF reader/writer (via vendored
  miniz), with zip-bomb / size caps on untrusted input.

## Architecture and contributing

- [`docs/algorithms.md`](docs/algorithms.md) — per-stage explanation
  of the repair pipeline.
- [`docs/internals.md`](docs/internals.md) — data structures, I/O
  formats, numerical conventions, build system.
- [`docs/README.md`](docs/README.md) — orientation for new contributors.

Pull requests welcome. Run the test suite locally before pushing:

```sh
cmake --build build --target tests_meshseal
./build/tests/tests_meshseal  --reporter=compact
# (Windows multi-config: ./build/tests/Release/tests_meshseal.exe)
```

## Vendored dependencies

All built from source via CMake FetchContent (no system packages
needed):

- [manifold](https://github.com/elalish/manifold) (Apache 2.0) — CSG
  boolean operations and level-set extraction.
- [miniz](https://github.com/richgel999/miniz) (MIT / public-domain)
  — ZIP archive reading for 3MF.
- [fast_float](https://github.com/fastfloat/fast_float) (Apache 2.0 /
  MIT / BSL) — locale-independent `from_chars<double>` polyfill.
- [Catch2](https://github.com/catchorg/Catch2) (BSL 1.0) — unit tests.

## Planned integrations

- **PrusaSlicer**: cross-platform alternative to the Windows-only
  repair backend currently in PrusaSlicer's GUI. Integration plan in
  the parent project's `docs/prusaslicer_integration_plan.md`.

## License

MIT — see [`LICENSE`](LICENSE).
