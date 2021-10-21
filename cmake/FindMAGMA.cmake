find_path(MAGMA_INCLUDE_DIR
  magma.h
  PATHS ${MAGA_ROOT} $ENV{MAGMA_ROOT}
  PATH_SUFFIXES include
)
mark_as_advanced(MAGMA_INCLUDE_DIR)

find_library(MAGMA_LIBRARIES
  NAMES magma
  PATHS ${MAGMA_ROOT} $ENV{MAGMA_ROOT}
  PATH_SUFFIXES lib
)
mark_as_advanced(MAGMA_LIBRARIES)

find_package_handle_standard_args(MAGMA DEFAULT_MSG MAGMA_LIBRARIES MAGMA_INCLUDE_DIR)

if(MAGMA_FOUND)
  find_package_message(MAGMA
    "Found MAGMA: ${MAGMA_LIBRARIES}"
    "[${MAGMA_LIBRARIES}][${MAGMA_INCLUDE_DIR}]"
  )
endif()

add_library(MAGMA INTERFACE)
target_include_directories(MAGMA INTERFACE ${MAGMA_INCLUDE_DIR})
target_link_libraries(MAGMA INTERFACE ${MAGMA_LIBRARIES})
