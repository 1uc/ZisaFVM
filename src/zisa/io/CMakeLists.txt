target_sources(zisa
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/backtrace.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/exec.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/file_manipulation.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/hdf5_serial_writer.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/hdf5_writer.cpp
)

