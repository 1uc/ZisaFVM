target_sources(unit_tests
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/barycentric.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/cartesian.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/cone.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/convergence_rates.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/edge.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/edge_rule.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/gauss_legendre.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/linear_interpolation.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/poly2d.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/quadrature.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/quasi_newton.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/reference_solution.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/subnormal.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/symmetric_coices.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/tetrahedron.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/triangle.cpp
)

