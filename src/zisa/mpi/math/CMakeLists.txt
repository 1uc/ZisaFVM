if(ZISA_HAS_MPI)

target_sources(zisa_mpi_obj
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/distributed_reference_solution.cpp
)

endif()
