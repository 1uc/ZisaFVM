if(ZISA_HAS_MPI)

endif()
if(ZISA_HAS_CUDA)

target_sources(zisa_cuda_obj
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/hello_world.cu
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/magma_context.cpp
)

endif()
