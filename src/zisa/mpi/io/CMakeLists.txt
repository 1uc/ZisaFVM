if(ZISA_HAS_MPI)

target_sources(zisa_mpi_obj
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/gathered_vis_info.cpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/gathered_visualization_factory.cpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/hdf5_unstructured_writer.cpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/mpi_progress_bar.cpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/parallel_load_snapshot.cpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/scattered_data_source_factory.cpp
)

endif()
