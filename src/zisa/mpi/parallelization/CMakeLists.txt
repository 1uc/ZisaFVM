if(ZISA_HAS_MPI)

target_sources(zisa_mpi
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/distributed_array_info.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/mpi_all_reduce.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/mpi_all_variables_gatherer.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/mpi_halo_exchange.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/mpi_single_node_array_gatherer.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/mpi_single_node_array_scatterer.cpp
)

endif()
