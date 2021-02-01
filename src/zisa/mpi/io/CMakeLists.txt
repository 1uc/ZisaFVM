if(ZISA_HAS_MPI)

target_sources(zisa_mpi
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/gathered_vis_info.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/gathered_visualization_factory.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/hdf5_unstructured_writer.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/mpi_progress_bar.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/parallel_load_snapshot.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/scattered_data_source_factory.cpp
)

endif()
add_subdirectory(cmake-build-debug)
