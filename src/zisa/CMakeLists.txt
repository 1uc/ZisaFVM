target_sources(zisa
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/datetime.cpp
)

add_subdirectory(cli)
add_subdirectory(flux)
add_subdirectory(grid)
add_subdirectory(io)
add_subdirectory(math)
add_subdirectory(model)
add_subdirectory(ode)
add_subdirectory(opengl)
add_subdirectory(reconstruction)
add_subdirectory(utils)
