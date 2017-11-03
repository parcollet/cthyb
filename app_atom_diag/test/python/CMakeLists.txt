# Copy h5 files to binary dir
file(GLOB all_h5_files *.h5)
file(COPY ${all_h5_files} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

# Get all .py files
file(GLOB all_tests *.py)

foreach(f ${all_tests})
  # t is the basename of the file
  get_filename_component(t ${f} NAME_WE)
  #
  add_test(NAME ${t} COMMAND python ${f})
  # Set the PythonPath : put the build dir first (in case there is an installed version). 
  set_property(TEST ${t} PROPERTY ENVIRONMENT PYTHONPATH=${CMAKE_BINARY_DIR}/python:$ENV{PYTHONPATH} )
endforeach()

