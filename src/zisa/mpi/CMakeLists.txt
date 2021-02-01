if(ZISA_HAS_MPI)

target_sources(zisa_mpi
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/mpi.cpp
)

endif()
add_subdirectory(io)
add_subdirectory(math)
add_subdirectory(parallelization)
