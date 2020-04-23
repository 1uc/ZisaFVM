if(ZISA_HAS_MPI)

target_sources(zisa_mpi_obj
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/hdf5_unstructured_writer.cpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/mpi_progress_bar.cpp
)

endif()
