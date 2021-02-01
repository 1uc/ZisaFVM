target_sources(zisa_generic_obj
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/barycentric.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/bounding_box.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/cartesian.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/cell.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/cell_factory.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/cone.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/denormalized_rule.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/edge.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/edge_rule.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/face.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/face_factory.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/geometric_spacing.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/linear_spacing.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/newton.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/permutation.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/poly2d.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/reference_solution.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/tetrahedral_rule.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/tetrahedron.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/triangle.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/triangular_rule.cpp
)

if(ZISA_HAS_MPI)

endif()
