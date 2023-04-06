add_executable(
    cg_hgp_hrr_path_gen
    cg_hgp_hrr_path_gen.cpp
)

target_link_libraries(
    cg_hgp_hrr_path_gen PRIVATE
    cg_base
)

add_executable(
    cg_hgp_hrr_code_gen
    cg_hgp_hrr_code_gen.cpp
)

target_link_libraries(
    cg_hgp_hrr_code_gen PRIVATE
    cg_base
)