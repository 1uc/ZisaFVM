target_sources(unit_tests
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/cweno_ao.cpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/lsq_solver.cpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/stencil.cpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/stencil_family.cpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/well_balanced_reconstruction.cpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/weno_ao.cpp
)

