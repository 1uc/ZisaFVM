target_sources(zisa_generic_obj
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/backtrace.cpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/colors.cpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/exec.cpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/file_manipulation.cpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/file_name_generator.cpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/gathered_visualization.cpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/hdf5_file.cpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/hdf5_serial_writer.cpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/hdf5_writer.cpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/no_visualization.cpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/progress_bar.cpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/scalar_plot.cpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/tri_plot.cpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/visualization.cpp
)

if(ZISA_HAS_MPI)
endif()
