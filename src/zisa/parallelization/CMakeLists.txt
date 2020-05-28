target_sources(zisa_generic_obj
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/all_reduce.cpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/all_variables_gatherer.cpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/all_variables_scatterer.cpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/distributed_grid.cpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/domain_decomposition.cpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/halo_info.cpp
)

if(ZISA_HAS_MPI)

endif()
