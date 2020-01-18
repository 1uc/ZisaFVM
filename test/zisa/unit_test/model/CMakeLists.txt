target_sources(unit_tests
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/all_variables.cpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/eos.cpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/janka.cpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/janka_eos.cpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/local_equilibrium.cpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/radial_poisson_solver.cpp
)

