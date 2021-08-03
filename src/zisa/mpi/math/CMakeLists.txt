if(ZISA_HAS_MPI)

target_sources(zisa_mpi_obj
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/distributed_reference_solution.cpp
)

endif()
if(ZISA_HAS_CUDA)

endif()
