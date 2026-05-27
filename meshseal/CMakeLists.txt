add_library(meshseal STATIC
    src/repair.cpp
    src/pipeline.cpp
    src/mesh.cpp
    src/stl_io.cpp
    src/3mf_io.cpp
    src/stages/weld.cpp
    src/stages/degenerate.cpp
    src/stages/orient.cpp
    src/stages/nm_vertex.cpp
    src/stages/nm_edge.cpp
    src/stages/non_manifold.cpp
    src/stages/holes.cpp
    src/stages/shells.cpp
    src/stages/intersections.cpp
    src/stages/thin_features.cpp
    src/stages/soup_reconstruct.cpp
    src/stages/fwn_levelset.cpp
    src/stages/sliver.cpp
    src/stages/nm_patch_remesh.cpp
    src/stages/bridge_loops.cpp
    src/stages/collapse_nm.cpp
    src/stages/tjunction.cpp
    src/stages/alpha_wrap.cpp
    src/stages/nm_local_repair.cpp
    src/stages/si_split.cpp
    src/stages/nm_carve_refill.cpp
    src/stages/coplanar_fan_drop.cpp
    src/stages/orient_wn.cpp
    src/stages/strip_doubled_membrane.cpp
    src/internal/winding_number.cpp
    src/internal/halfedge.cpp
    src/internal/bvh.cpp
    src/internal/diagnostics.cpp
    src/internal/component_classifier.cpp
    src/internal/voxel_grid.cpp
    src/internal/arrangement.cpp
    src/internal/si_count.cpp
)

target_include_directories(meshseal
    PUBLIC  include
    PRIVATE src
)

target_link_libraries(meshseal
    PUBLIC
        manifold
        miniz
    PRIVATE
        fast_float
)

target_compile_options(meshseal PRIVATE
    $<$<CXX_COMPILER_ID:MSVC>:/W4>
    $<$<OR:$<CXX_COMPILER_ID:GNU>,$<CXX_COMPILER_ID:Clang>>:-Wall>
    $<$<OR:$<CXX_COMPILER_ID:GNU>,$<CXX_COMPILER_ID:Clang>>:-Wextra>
)
