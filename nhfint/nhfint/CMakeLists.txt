add_subdirectory(src/boysfun)       # boys function
add_subdirectory(src/ryspoly)       # roots and weights of rys polynomial
add_subdirectory(src/autocode)      # recurrence relations for md, dkr and hgp

add_library(
    nhfint
    src/angmom_cart.cpp
    src/basisfile.cpp
    src/utility.cpp
    src/basis.cpp
    src/int1e_class_cart_d0.cpp
    src/int1e_class_cart_d1.cpp
    src/int1e_class_cart_d2.cpp
    src/int2e_class_cart_d0.cpp
    src/int2e_class_cart_d1.cpp
    src/int2e_class_cart_d2.cpp
    src/int1e_tensor_d0.cpp
    src/int1e_tensor_d1.cpp
    src/int1e_tensor_d2.cpp
    src/int2e_tensor_d0.cpp
    src/int2e_tensor_d1.cpp
    src/int2e_tensor_d2.cpp
    src/calc_ericlass_d0_cart_dkr.cpp
    src/calc_ericlass_d0_cart_hgp.cpp
    src/calc_ericlass_d0_cart_tho.cpp
)

target_link_libraries(
    nhfint PUBLIC
    boysfun                 # boys function
    ryspoly                 # roots and weights of rys polynomial
    nhfmath
    nhfstr
)

target_include_directories(
    nhfint PUBLIC
    include
)


# test
add_subdirectory(test)