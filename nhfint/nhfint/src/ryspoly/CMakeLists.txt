# generate library
add_library(
    ryspoly
    src/roots_weights.cpp
)

target_include_directories(
    ryspoly PUBLIC
    include
)

