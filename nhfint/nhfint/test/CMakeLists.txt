# testing of basisfile
add_custom_target(
  copy_test_data_basisfile
  ${CMAKE_COMMAND} -E create_symlink
  ${CMAKE_CURRENT_SOURCE_DIR}/test_data_basisfile
  ${CMAKE_CURRENT_BINARY_DIR}/test_data_basisfile
)

add_executable(
    test_basisfile
    test_basisfile.cpp
)

target_link_libraries(
    test_basisfile PRIVATE
    nhfint
    gtest
    gtest_main
)

add_dependencies(
    test_basisfile
    copy_test_data_basisfile
)


# testing angmom_cart
add_executable(
    test_angmom_cart
    test_angmom_cart.cpp
)

target_link_libraries(
    test_angmom_cart
    nhfint
    gtest
    gtest_main
)


# testing utility
add_executable(
    test_utility
    test_utility.cpp
)

target_link_libraries(
    test_utility
    nhfint
    gtest
    gtest_main
)


# # testing calc_ericlass_d0_dkr
# add_executable(
#     test_calc_ericlass_d0_cart_dkr
#     test_calc_ericlass_d0_cart_dkr.cpp
# )

# target_link_libraries(
#     test_calc_ericlass_d0_cart_dkr
#     nhfint
#     gtest
#     gtest_main
# )


# testing calc_ericlass_d0_dkr
add_executable(
    test_calc_ericlass_d0_cart_tho
    test_calc_ericlass_d0_cart_tho.cpp
)

target_link_libraries(
    test_calc_ericlass_d0_cart_tho
    nhfint
    gtest
    gtest_main
)