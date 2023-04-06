add_library(
    nhfstr
    src/nhfstr.cpp
)

target_include_directories(
    nhfstr PUBLIC
    include
)


add_subdirectory(test)