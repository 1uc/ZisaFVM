target_sources(zisa_generic
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/backtrace.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/colors.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/data_source.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/exec.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/file_name_generator.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/gathered_visualization.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/load_snapshot.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/no_visualization.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/progress_bar.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/scalar_plot.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/scattered_data_source.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/tri_plot.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/visualization.cpp
)

if(ZISA_HAS_MPI)

endif()
