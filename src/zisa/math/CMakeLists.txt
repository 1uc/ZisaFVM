target_sources(zisa_obj
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/barycentric.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/basic_functions.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/bounding_box.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/cone.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/edge.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/edge_rule.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/newton.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/poly2d.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/reference_solution.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/triangle.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/triangular_rule.cpp
)

