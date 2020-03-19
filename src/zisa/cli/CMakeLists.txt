target_sources(zisa_generic_obj
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/input_parameters.cpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/parse_command_line.cpp
)

if(ZISA_HAS_MPI)
endif()
