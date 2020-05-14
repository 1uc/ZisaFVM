add_subdirectory(zisa)

target_sources(zisa
  PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/zisa.cpp
)

target_sources(locate-point
  PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/locate_point.cpp
)

target_sources(domain-decomposition
  PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/domain_decomposition.cpp
)

target_sources(opengl-demo
  PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/opengl_demo.cpp
)

target_sources(something-demo
  PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/something_demo.cpp
)

target_sources(space-filling-curve-demo
  PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/space_filling_curve_demo.cpp
)
