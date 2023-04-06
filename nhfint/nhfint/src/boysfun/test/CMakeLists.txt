# testing
add_custom_target(
  copy_test_data
  ${CMAKE_COMMAND} -E create_symlink
  ${CMAKE_CURRENT_SOURCE_DIR}/test_data
  ${CMAKE_CURRENT_BINARY_DIR}/test_data
)

add_executable(
    test_boysfun
    test_boysfun.cpp
)

target_link_libraries(
    test_boysfun PRIVATE
    boysfun
    gtest
    gtest_main
)

add_dependencies(
    test_boysfun
    copy_test_data
)


# benchmark
add_executable(
    bench_boysfun
    bench_boysfun.cpp
)

target_link_libraries(
    bench_boysfun PRIVATE
    boysfun
    benchmark
    benchmark_main
    pthread
)

add_dependencies(
    bench_boysfun
    copy_test_data
)
