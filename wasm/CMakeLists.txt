# WebAssembly build target for meshseal.
#
# Only added when MESHSEAL_BUILD_WASM=ON, which is only meaningful under the
# Emscripten toolchain (configure with `emcmake cmake … -DMESHSEAL_BUILD_WASM=ON
# -DMESHSEAL_BUILD_CLI=OFF -DMESHSEAL_BUILD_TESTS=OFF`). Produces meshseal.js +
# meshseal.wasm: a MODULARIZE'd ES-factory loadable from a Web Worker or Node.

if(NOT EMSCRIPTEN)
    message(FATAL_ERROR
        "MESHSEAL_BUILD_WASM requires the Emscripten toolchain. "
        "Configure with: emcmake cmake -S . -B build-wasm -DMESHSEAL_BUILD_WASM=ON")
endif()

add_executable(meshseal_wasm meshseal_wasm.cpp)

target_link_libraries(meshseal_wasm PRIVATE meshseal)

# The shim uses the private buffer-I/O header (src/stl_io.h). The public
# meshseal.h comes in transitively via the meshseal target's PUBLIC include.
target_include_directories(meshseal_wasm PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/../meshseal/src)

# Output as meshseal.js (+ meshseal.wasm beside it).
set_target_properties(meshseal_wasm PROPERTIES
    OUTPUT_NAME "meshseal"
    SUFFIX ".js")

# embind for the JS<->C++ binding; MODULARIZE so the loader is a factory
# function (createMeshseal()) rather than polluting global scope — required
# for clean Web Worker + Node use. ALLOW_MEMORY_GROWTH because mesh size is
# unknown up front. ENVIRONMENT spans node/web/worker so the same artifact
# runs in the spike harness (Node) and later in the browser worker.
target_link_options(meshseal_wasm PRIVATE
    "-lembind"
    "-O2"
    "-sMODULARIZE=1"
    "-sEXPORT_NAME=createMeshseal"
    "-sALLOW_MEMORY_GROWTH=1"
    "-sENVIRONMENT=node,web,worker"
    "-sEXIT_RUNTIME=0"
    # Surface C++ exceptions to embind try/catch (the shim catches
    # std::exception and reports {ok:false,error}). Without this, a throw
    # would abort the module instead of returning an error object.
    "-fexceptions")

# C++ exceptions must be enabled at compile time too for -fexceptions to work.
target_compile_options(meshseal_wasm PRIVATE "-fexceptions")
