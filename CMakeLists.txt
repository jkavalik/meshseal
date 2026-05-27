cmake_minimum_required(VERSION 3.18)
project(meshseal VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Build the FetchContent deps statically so the CLI / test binaries have
# no runtime DLL dependencies. Without this:
#   - Windows: catch_discover_tests fails with 0xc0000135 because
#     tests_meshseal.exe can't find manifold.dll / miniz.dll at build time.
#   - Linux/macOS: standalone meshseal_cli requires libmanifold.so / .dylib
#     on LD_LIBRARY_PATH / DYLD_LIBRARY_PATH at runtime.
# Static deps make every binary self-contained.
set(BUILD_SHARED_LIBS OFF)

# miniz 3.0.2 declares cmake_minimum_required(VERSION 3.0) at the top of
# its own CMakeLists. CMake 4+ refuses any minimum below 3.5 unless this
# policy-version cap is raised. Affects macOS runners that have shipped
# CMake 4.x already and any contributor with a recent CMake install.
set(CMAKE_POLICY_VERSION_MINIMUM 3.5)

include(CTest)

include(FetchContent)


FetchContent_Declare(
    Catch2
    GIT_REPOSITORY https://github.com/catchorg/Catch2.git
    GIT_TAG        v3.7.1
)

FetchContent_Declare(
    manifold
    GIT_REPOSITORY https://github.com/elalish/manifold.git
    GIT_TAG        v3.4.1
)

FetchContent_Declare(
    fast_float
    GIT_REPOSITORY https://github.com/fastfloat/fast_float.git
    GIT_TAG        v6.1.6
    # Polyfill for std::from_chars<double>. The C++17 standard requires
    # the floating-point overload, but real-world toolchains lag:
    #   - GCC 9, 10: no double overload (compile error)
    #   - libc++ in Apple Clang up through Xcode 16: overload explicitly
    #     `= delete`d (compile error on macos-latest CI)
    #   - libc++ < 14 on Linux: same issue
    # fast_float is a header-only, MIT/Apache/BSL triple-licensed library
    # that exposes fast_float::from_chars with the same API. The integer
    # overload of std::from_chars is universally supported (we keep it).
)

FetchContent_Declare(
    miniz
    GIT_REPOSITORY https://github.com/richgel999/miniz.git
    GIT_TAG        3.0.2
    # NOTE on version reporting: the fetched miniz advertises three
    # mutually-inconsistent version strings — `miniz.c 3.0.0` in the source
    # header comment, `GIT_TAG 3.0.2` here, and `MZ_VERSION "11.0.2"` macro.
    # This is upstream's idiosyncrasy (MZ_VERSION is legacy internal
    # numbering that hasn't been bumped in years); the GIT_TAG is what
    # actually pins the build. If supply-chain robustness becomes a
    # concern, lock by commit SHA instead of tag string — `git ls-remote
    # https://github.com/richgel999/miniz.git refs/tags/3.0.2` to obtain
    # the SHA, then GIT_TAG <sha>.
)

set(MANIFOLD_TEST             OFF CACHE BOOL "" FORCE)
set(MANIFOLD_PYBIND           OFF CACHE BOOL "" FORCE)
set(MANIFOLD_CBIND            OFF CACHE BOOL "" FORCE)
set(MANIFOLD_FUZZ             OFF CACHE BOOL "" FORCE)
set(MANIFOLD_EXPORT           OFF CACHE BOOL "" FORCE)
set(MANIFOLD_BUILD_BENCHMARKS OFF CACHE BOOL "" FORCE)
set(MANIFOLD_BUILD_EXAMPLES   OFF CACHE BOOL "" FORCE)

FetchContent_MakeAvailable(Catch2 manifold miniz fast_float)

add_subdirectory(meshseal)
add_subdirectory(cli)
add_subdirectory(tests)
