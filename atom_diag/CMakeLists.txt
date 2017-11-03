file(GLOB HEADERS *.hpp)
file(GLOB_RECURSE SOURCES *.cpp)

add_library(atom_diag_c ${SOURCES})
target_link_libraries(atom_diag_c PUBLIC triqs)
#triqs_set_rpath_for_target(atom_diag_c)

# Install
install(TARGETS atom_diag_c DESTINATION lib)
install(FILES ${HEADERS} DESTINATION include/triqs/atom_diag)

