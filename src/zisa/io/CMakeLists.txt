target_sources(zisa
  PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/backtrace.cpp
  PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/exec.cpp
  PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/file_manipulation.cpp
  PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/hdf5_writer.cpp
  PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/hdf5_serial_writer.cpp
)
