# generate library
add_library(
    boysfun
    src/boysfun.cpp
)

target_include_directories(
    boysfun PUBLIC
    include
)


# testing and benchmark
add_subdirectory(test)