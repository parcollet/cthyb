# Build the python module
add_cpp2py_module(atom_diag)
target_link_libraries(atom_diag atom_diag_c)

# We need to include the convertes.hxx files
target_include_directories(atom_diag PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

# Add version and __init__ file
#configure_file(version.py.in version.py)

set(PYTHON_SOURCES
 __init__.py
 #${CMAKE_CURRENT_BINARY_DIR}/version.py
 )

foreach(f ${PYTHON_SOURCES})
 configure_file(${f} ${f} COPYONLY)
endforeach()

# Install python module to proper location
set(PYTHON_LIB_DEST ${CPP2PY_PYTHON_LIB_DEST_ROOT}/triqs/atom_diag)
install(TARGETS atom_diag DESTINATION ${PYTHON_LIB_DEST})
install(FILES ${PYTHON_SOURCES} DESTINATION ${PYTHON_LIB_DEST})

