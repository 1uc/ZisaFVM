target_sources(zisa_generic_obj
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/all_variables.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/distributed_cfl_condition.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/euler.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/euler_factory.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/euler_variables.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/gravity.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/grid_variables.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/helmholtz_eos.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/ideal_gas_eos.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/instantaneous_physics.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/janka_eos.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/load_full_state.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/models.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/polytrope.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/radial_poisson_solver.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/save_full_state.cpp
)

if(ZISA_HAS_MPI)

endif()
