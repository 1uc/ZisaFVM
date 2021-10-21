if(ZISA_HAS_MPI)

endif()
if(ZISA_HAS_CUDA)

target_sources(zisa_cuda_obj
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/cholesky_decomposition.cpp
)

endif()
