set(LIBRARY_SOURCES
    solver_core.cpp
    atom_diag.cpp
    atom_diag_functions.cpp
    atom_diag_worker.cpp
    impurity_trace.cpp
    moves/insert.cpp
    moves/remove.cpp
    moves/shift.cpp
    moves/global.cpp
    moves/double_insert.cpp
    moves/double_remove.cpp
    measures/density_matrix.cpp
    measures/G_tau.cpp
    measures/G_l.cpp
    )

#FIXME : for cmake > 3.1, use target_sources below
if(MeasureG2)
  list(APPEND LIBRARY_SOURCES
       measures/G2_tau.cpp
       measures/G2_iw.cpp
       measures/G2_iwll.cpp
       )
endif()

add_library(cthyb_c ${LIBRARY_SOURCES})
target_link_libraries(cthyb_c PUBLIC triqs)

# FIXME : should go ?
option(EXT_DEBUG "Enable extended debugging output [developers only]" OFF)

if(EXT_DEBUG)
 target_compile_options(cthyb_c PRIVATE -DEXT_DEBUG)
endif()

# FIXME : should go
option(SAVE_CONFIGS "Save visited configurations to configs.h5 [developers only]" OFF)
if(SAVE_CONFIGS)
 set(NUM_CONFIGS_TO_SAVE 50000 CACHE STRING "Number of visited configurations to save [developers only]")
 target_compile_options(cthyb_c PRIVATE -DSAVE_CONFIGS -DNUM_CONFIGS_TO_SAVE=${NUM_CONFIGS_TO_SAVE})
endif()

# 2-particle GF measurement requires NFFT
if(MeasureG2)
 target_link_libraries(cthyb_c PRIVATE nfft)
 #target_sources(cthyb_c PRIVATE measures/G2_iw.cpp measures/G2_iwll.cpp)
endif()

# Install
install(TARGETS cthyb_c DESTINATION lib)

#triqs_set_rpath_for_target(cthyb_c)

# FIXME : WHY DO WE INSTALL HEADERS HERE ??
set(HEADERS_TO_INSTALL
    array_suppl.hpp # FIXME
    atom_diag.hpp
    atom_diag_functions.hpp
    solve_parameters.hpp
    solver_core.hpp
    config.hpp
)
install(FILES ${HEADERS_TO_INSTALL} DESTINATION include/triqs_cthyb)
