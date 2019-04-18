target_sources(unit_tests
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/all_variables.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/eos.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/janka.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/janka_eos.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/local_equilibrium.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/radial_poisson_solver.cpp
)

