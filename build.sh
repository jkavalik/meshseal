#!/usr/bin/env bash
# build.sh — build the meshseal library, CLI, and test binary.
#
# Usage (from the repository root or from library/):
#   ./library/build.sh           # incremental build, all targets
#   ./library/build.sh --clean   # wipe build tree, reconfigure, then build
#   ./library/build.sh --test    # build + run tests via CTest
#   ./library/build.sh --clean --test

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SOURCE_DIR="$SCRIPT_DIR"
BUILD_DIR="$SCRIPT_DIR/build"

DO_CLEAN=0
DO_TEST=0
for arg in "$@"; do
    case "$arg" in
        --clean|-c) DO_CLEAN=1 ;;
        --test|-t)  DO_TEST=1  ;;
        *) echo "Unknown argument: $arg"; exit 1 ;;
    esac
done

if [[ $DO_CLEAN -eq 1 && -d "$BUILD_DIR" ]]; then
    echo "-- Removing build tree: $BUILD_DIR"
    rm -rf "$BUILD_DIR"
fi

if [[ ! -d "$BUILD_DIR" ]]; then
    echo "-- Configuring CMake..."
    cmake -S "$SOURCE_DIR" -B "$BUILD_DIR" -DCMAKE_BUILD_TYPE=Release
fi

echo "-- Building all targets (Release)..."
cmake --build "$BUILD_DIR" --config Release --target meshseal_cli tests_meshseal

echo "-- Build succeeded."
echo "   CLI:   $BUILD_DIR/cli/Release/meshseal_cli.exe"
echo "   Tests: $BUILD_DIR/tests/Release/tests_meshseal.exe"

if [[ $DO_TEST -eq 1 ]]; then
    echo ""
    echo "-- Running CTest..."
    cd "$BUILD_DIR"
    ctest -C Release --output-on-failure
fi
