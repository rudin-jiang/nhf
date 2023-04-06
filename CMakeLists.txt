cmake_minimum_required(VERSION 3.10)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

project(nhf C CXX)

if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

add_subdirectory(extern/googletest)
add_subdirectory(extern/benchmark)

add_subdirectory(nhfstr)
add_subdirectory(nhfdata)
add_subdirectory(nhfmath)
add_subdirectory(nhfint)
add_subdirectory(nhf)