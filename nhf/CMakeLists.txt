add_executable(
    nhf
    src/rhf.cpp
    src/nhf_main.cpp
)

target_link_libraries(
    nhf PRIVATE
    nhfmath
    nhfstr
    nhfint
)

target_include_directories(
    nhf PRIVATE
    include
)






add_subdirectory(test)