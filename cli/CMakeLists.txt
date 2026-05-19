add_executable(meshseal_cli main.cpp)

target_link_libraries(meshseal_cli PRIVATE meshseal)

install(TARGETS meshseal_cli
    RUNTIME DESTINATION bin
)
