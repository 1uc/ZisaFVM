target_sources(zisa_generic
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/all_reduce.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/all_variables_gatherer.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/all_variables_scatterer.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/distributed_grid.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/domain_decomposition.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/halo_exchange.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/halo_info.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/local_grid.cpp
)

if(ZISA_HAS_MPI)

endif()
