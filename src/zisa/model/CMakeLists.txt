target_sources(zisa_obj
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/all_variables.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/euler_factory.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/euler_variables.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/gravity.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/grid_variables.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/ideal_gas_eos.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/instantaneous_physics.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/janka_eos.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/polytrope.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/radial_poisson_solver.cpp
)

