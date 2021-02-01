target_sources(zisa_generic_obj
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/load_shaders.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/tri_plot.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/window.cpp
)

if(ZISA_HAS_MPI)

endif()
