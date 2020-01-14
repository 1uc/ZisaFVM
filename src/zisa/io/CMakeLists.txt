target_sources(zisa_obj
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/backtrace.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/colors.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/exec.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/file_manipulation.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/file_name_generator.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/gathered_visualization.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/hdf5_file.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/hdf5_serial_writer.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/mpi_progress_bar.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/no_visualization.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/progress_bar.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/scalar_plot.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/tri_plot.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/visualization.cpp
)

