target_sources(unit_tests
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/magma_batched.cpp
)

add_subdirectory(reconstruction)
