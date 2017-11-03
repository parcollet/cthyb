# Common name
set (APP atom_diag)

# Version number 
set (ATOM_DIAG_VERSION "1.5")

# Start configuration
cmake_minimum_required(VERSION 2.8)
project(triqs_${APP} CXX)

# Default to Release build type
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Type of build" FORCE)
endif()
message( STATUS "-------- BUILD-TYPE: ${CMAKE_BUILD_TYPE} -------------")

# Use shared libraries
set(BUILD_SHARED_LIBS ON)

# Load TRIQS and Cpp2Py
find_package(TRIQS 1.5 REQUIRED)
find_package(Cpp2Py REQUIRED)

# Default Install directory to TRIQS_ROOT
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT OR (NOT IS_ABSOLUTE ${CMAKE_INSTALL_PREFIX}))
  message(STATUS " No install prefix given (or invalid). Defaulting to TRIQS_ROOT")
  set(CMAKE_INSTALL_PREFIX ${TRIQS_ROOT} CACHE PATH "default install path" FORCE)
endif()
message(STATUS "-------- CMAKE_INSTALL_PREFIX: ${CMAKE_INSTALL_PREFIX} -------------")

# Macro defined in TRIQS which picks the hash of repo.
triqs_get_git_hash_of_source_dir(ATOM_DIAG_GIT_HASH)

# Options for compilation that apply to everythings (C++ solver, Python interface, tests).
add_compile_options( -std=c++14 -DHAS_OPTIONAL_HEADER 
                     -DATOM_DIAG_GIT_HASH=${ATOM_DIAG_GIT_HASH} -DTRIQS_GIT_HASH=${TRIQS_GIT_HASH}
                      $<$<CONFIG:Release>:-march=native>
                      )

# Common to all targets below (C++ lib, python, tests).
include_directories(${CMAKE_SOURCE_DIR})

# To be fixed!!
triqs_set_rpath_for_all_targets()

add_subdirectory(${APP})

# Tests
enable_testing()
add_subdirectory(test)

if (${TRIQS_WITH_PYTHON_SUPPORT})

 # Python interface
 add_subdirectory(python/${APP})
 
 # Build the documentation
 option(Build_Documentation "Build documentation" OFF)
 if (${Build_Documentation})
    add_subdirectory(doc)
 endif()
 
endif()

