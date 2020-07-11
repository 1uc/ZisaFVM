target_sources(zisa_generic_obj
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/all_variables.cpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/distributed_cfl_condition.cpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/euler.cpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/euler_factory.cpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/euler_variables.cpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/gravity.cpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/grid_variables.cpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/helmholtz_eos.cpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/ideal_gas_eos.cpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/instantaneous_physics.cpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/janka_eos.cpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/models.cpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/polytrope.cpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/radial_poisson_solver.cpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/save_full_state.cpp
)

if(ZISA_HAS_MPI)

endif()
