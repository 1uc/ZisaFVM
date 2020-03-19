target_sources(zisa_generic_obj
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/cweno_ao.cpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/hybrid_weno.cpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/hybrid_weno_params.cpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/lsq_solver.cpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/lsq_solver_family.cpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/stencil.cpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/stencil_bias.cpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/stencil_family.cpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/stencil_family_params.cpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/stencil_params.cpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/weno_ao.cpp
)

if(ZISA_HAS_MPI)
endif()
