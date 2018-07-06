add_subdirectory(zisa)

target_sources(run
  PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/zisa.cpp
)
