include("cmake/GetCompilerVersion.cmake")
get_filename_component(C_COMPILER_NAME ${CMAKE_C_COMPILER} NAME)

set(METIS_VER 5.1.0)
set(METIS_DIR ${PROJECT_SOURCE_DIR}/third_party/metis-${METIS_VER}/${C_COMPILER_NAME}/${C_COMPILER_VERSION})

add_library(metis INTERFACE)
target_include_directories(metis INTERFACE ${METIS_DIR}/include)
target_link_directories(metis INTERFACE ${METIS_DIR}/lib)
target_link_libraries(metis INTERFACE -lmetis)