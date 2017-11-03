# Version number 
set (CTHYB_VERSION "1.5")

# start configuration
cmake_minimum_required(VERSION 2.8)
project(cthyb CXX)

# Default to Release build type
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Type of build" FORCE)
endif()
message( STATUS "-------- BUILD-TYPE: ${CMAKE_BUILD_TYPE} -------------")

# Options for the compilation ...
option(Hybridisation_is_complex "If ON, the hybridization Delta(tau) is complex" OFF)
option(Local_hamiltonian_is_complex "If ON, the H_loc is complex" OFF)
option(MeasureG2 "Measure the two particle object (requires the NFFT library)" ON)

# check that options are compatible
if(Hybridisation_is_complex AND NOT Local_hamiltonian_is_complex)
 message(FATAL_ERROR "Combination of options Hybridisation_is_complex=ON "
                     "and Local_hamiltonian_is_complex=OFF is not supported.")
endif()

# Use shared libraries
set(BUILD_SHARED_LIBS ON)

# Load TRIQS and Cpp2Py
find_package(TRIQS 1.5 EXACT REQUIRED)
find_package(Cpp2Py REQUIRED)

# We will need the FindNNFT of triqs, so we adjust CMake Module Path to find it
list(APPEND CMAKE_MODULE_PATH ${TRIQS_ROOT}/share/cmake)

# 2-particle GF measurement requires NFFT
if(MeasureG2)
 find_package(NFFT REQUIRED)
endif()

# Default Install directory to TRIQS_ROOT if not given. Checks an absolute name is given.
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT OR (NOT IS_ABSOLUTE ${CMAKE_INSTALL_PREFIX}))
  message(STATUS " No install prefix given (or invalid). Defaulting to TRIQS_ROOT")
  set(CMAKE_INSTALL_PREFIX ${TRIQS_ROOT} CACHE PATH "default install path" FORCE)
endif()
message(STATUS "-------- CMAKE_INSTALL_PREFIX: ${CMAKE_INSTALL_PREFIX} -------------")

# Macro defined in TRIQS which picks the hash of repo.
triqs_get_git_hash_of_source_dir(CTHYB_GIT_HASH)
message(STATUS "Git hash: ${CTHYB_GIT_HASH}")

# Options for compilation that apply to everythings (C++ solver, Python interface, tests).
add_compile_options(-std=c++14 
                    -DCTHYB_GIT_HASH=${CTHYB_GIT_HASH} -DTRIQS_GIT_HASH=${TRIQS_GIT_HASH}
                     $<$<CONFIG:Release>:-march=native>
                     $<$<CONFIG:Debug>:-DDEBUG_CTHYB>
                     $<$<BOOL:${Hybridisation_is_complex}>:-DHYBRIDISATION_IS_COMPLEX>
                     $<$<BOOL:${Local_hamiltonian_is_complex}>:-DLOCAL_HAMILTONIAN_IS_COMPLEX>
                     $<$<BOOL:${MeasureG2}>:-DCTHYB_G2_NFFT>
                     )
# Common to all targets (ctint, python, tests).
include_directories(${CMAKE_SOURCE_DIR} ${CMAKE_SOURCE_DIR}/app_atom_diag)

# To be fixed!!
triqs_set_rpath_for_all_targets()

add_subdirectory(cthyb)
add_subdirectory(app_atom_diag/atom_diag)

# Tests
option(Tests "Enable Tests" ON)
if (${Tests})
 enable_testing()
 add_subdirectory(app_atom_diag/test)
 add_subdirectory(test)
endif()

if (${TRIQS_WITH_PYTHON_SUPPORT})

 # Python interface
 add_subdirectory(app_atom_diag/python/atom_diag)
 add_subdirectory(python/cthyb)

 # For cpp2py generation to see atom_diag in ./python
 execute_process(COMMAND ln -s ../app_atom_diag/python/atom_diag  python/atom_diag) 

 # Build the documentation
 option(Build_Documentation "Build documentation" OFF)
 if (${Build_Documentation})
   add_subdirectory(app_atom_diag/doc)
   add_subdirectory(doc)
 endif()
 
endif()




