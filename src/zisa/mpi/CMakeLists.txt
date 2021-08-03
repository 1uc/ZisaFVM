if(ZISA_HAS_MPI)

endif()
if(ZISA_HAS_CUDA)

endif()
add_subdirectory(io)
add_subdirectory(math)
add_subdirectory(parallelization)
