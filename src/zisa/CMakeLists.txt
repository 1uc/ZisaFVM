target_sources(zisa_obj
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/datetime.cpp
)

add_subdirectory(boundary)
add_subdirectory(cli)
add_subdirectory(core)
add_subdirectory(experiments)
add_subdirectory(flux)
add_subdirectory(grid)
add_subdirectory(io)
add_subdirectory(loops)
add_subdirectory(math)
add_subdirectory(model)
add_subdirectory(ode)
add_subdirectory(opengl)
add_subdirectory(reconstruction)
add_subdirectory(utils)
