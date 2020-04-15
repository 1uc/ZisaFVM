target_sources(zisa_generic_obj
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/all_reduce.cpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/all_variables_gatherer.cpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/distributed_grid.cpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/domain_decomposition.cpp
)

if(ZISA_HAS_MPI)
target_sources(zisa_mpi_obj
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/mpi.cpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/mpi_all_reduce.cpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/mpi_all_variables_gatherer.cpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/mpi_halo_exchange.cpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/mpi_single_node_array_gatherer.cpp
)

endif()
