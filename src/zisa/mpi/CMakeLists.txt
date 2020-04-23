if(ZISA_HAS_MPI)

target_sources(zisa_mpi_obj
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/mpi.cpp
)

endif()
add_subdirectory(io)
add_subdirectory(parallelization)
