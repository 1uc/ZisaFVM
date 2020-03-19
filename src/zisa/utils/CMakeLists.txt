target_sources(zisa_generic_obj
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/indent_block.cpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/logging.cpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/parse_duration.cpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/replace_all.cpp
)

if(ZISA_HAS_MPI)
endif()
