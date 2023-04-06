find_package(Eigen3 REQUIRED)

add_library(
    nhfmath
    src/vec3d.cpp
    src/matrix.cpp
    src/mathfun.cpp
)

target_include_directories(
    nhfmath PUBLIC
    include
    ${EIGEN3_INCLUDE_DIR}
)


# testing
add_subdirectory(test)