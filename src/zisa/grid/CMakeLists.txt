target_sources(zisa_generic_obj
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/cell.cpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/gmsh_reader.cpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/grid.cpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/point_locator.cpp
)

if(ZISA_HAS_MPI)

endif()
