target_sources(zisa_generic_obj
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/barycentric.cpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/basic_functions.cpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/bounding_box.cpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/cell.cpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/cell_factory.cpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/cone.cpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/denormalized_rule.cpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/edge.cpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/edge_rule.cpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/face.cpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/face_factory.cpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/geometric_spacing.cpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/linear_spacing.cpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/newton.cpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/permutation.cpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/poly2d.cpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/reference_solution.cpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/tetrahedral_rule.cpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/tetrahedron.cpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/triangle.cpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/triangular_rule.cpp
)

if(ZISA_HAS_MPI)
endif()
