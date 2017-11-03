# Copy h5 files to binary dir
file(GLOB all_h5_files *.h5)
file(COPY ${all_h5_files} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

# Get all .cpp files
file(GLOB all_tests *.cpp)

foreach(f ${all_tests})
  get_filename_component (t ${f} NAME_WE)
  add_executable(${t} ${f})
  target_link_libraries(${t} atom_diag_c gtest)
  add_test(${t} ${CMAKE_CURRENT_BINARY_DIR}/${t})
endforeach()
