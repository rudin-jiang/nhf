add_executable(
    cg_hgp_vrr_path_gen
    cg_hgp_vrr_path_gen.cpp
)

target_link_libraries(
    cg_hgp_vrr_path_gen PRIVATE
    cg_base
)