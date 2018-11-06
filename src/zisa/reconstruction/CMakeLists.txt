target_sources(zisa
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/cweno_ao.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/hybrid_weno.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/hybrid_weno_params.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/lsq_solver.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/lsq_solver_family.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/stencil.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/stencil_bias.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/stencil_family.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/stencil_family_params.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/stencil_params.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/weno_ao.cpp
)

