add_library(
    cg_base
    src/ijk.cpp
    src/ijkpair.cpp
    src/ijkpairm.cpp
    src/iterm.cpp
    src/utility.cpp
)

target_include_directories(
    cg_base PUBLIC
    include
)


# test
add_subdirectory(
    test
)