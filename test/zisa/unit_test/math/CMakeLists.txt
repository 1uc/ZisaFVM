target_sources(unit_tests
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/barycentric.cpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/cartesian.cpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/cone.cpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/convergence_rates.cpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/edge.cpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/edge_rule.cpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/face.cpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/gauss_legendre.cpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/linear_interpolation.cpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/permutation.cpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/poly2d.cpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/quadrature.cpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/quasi_newton.cpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/reference_solution.cpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/subnormal.cpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/symmetric_coices.cpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/tetrahedral_rule.cpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/tetrahedron.cpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/triangle.cpp
)

