target_sources(unit_tests
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/cweno_ao.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/lsq_solver.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/stencil.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/stencil_family.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/well_balanced_reconstruction.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/weno_ao.cpp
)

