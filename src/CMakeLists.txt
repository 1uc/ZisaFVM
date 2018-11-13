add_subdirectory(zisa)

target_sources(run
  PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/zisa.cpp
)

target_sources(opengl-demo
  PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/opengl_demo.cpp
)
