cmake_minimum_required(VERSION 3.21)
project(meshseal VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

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

FetchContent_MakeAvailable(Catch2 manifold miniz)

add_subdirectory(meshseal)
add_subdirectory(cli)
add_subdirectory(tests)
