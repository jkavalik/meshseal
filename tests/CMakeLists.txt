include(Catch)

add_executable(tests_meshseal
    test_stl_io.cpp
    test_3mf_io.cpp
    test_weld.cpp
    test_e2e.cpp
    test_degenerate.cpp
    test_orient.cpp
    test_nm_vertex.cpp
    test_nm_edge.cpp
    test_holes.cpp
    test_shells.cpp
    test_intersections.cpp
    test_thin_features.cpp
    test_soup_reconstruct.cpp
    test_pipeline_architecture.cpp
    test_alpha_wrap.cpp
    test_progress_callback.cpp
)

target_link_libraries(tests_meshseal
    PRIVATE
        Catch2::Catch2WithMain
        meshseal
)

catch_discover_tests(tests_meshseal)
