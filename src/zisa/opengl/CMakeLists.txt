target_sources(zisa_generic_obj
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/load_shaders.cpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/tri_plot.cpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/window.cpp
)

if(ZISA_HAS_MPI)
endif()
