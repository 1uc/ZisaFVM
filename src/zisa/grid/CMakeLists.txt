target_sources(zisa_generic_obj
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/cell.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/gmsh_reader.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/grid.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/point_locator.cpp
)

if(ZISA_HAS_MPI)

endif()
if(ZISA_HAS_CUDA)

endif()
