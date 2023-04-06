# test ijk
add_executable(
    test_ijk
    test_ijk.cpp
)

target_link_libraries(
    test_ijk PRIVATE
    cg_base
    gtest
    gtest_main
)