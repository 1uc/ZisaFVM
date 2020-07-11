include("cmake/GetCompilerVersion.cmake")
get_filename_component(C_COMPILER_NAME ${CMAKE_C_COMPILER} NAME)

set(HELMHOLTZ_EOS_DIR ${PROJECT_SOURCE_DIR}/third_party/helmholtz_eos/${C_COMPILER_NAME}/${C_COMPILER_VERSION})

add_library(helmholtz_eos INTERFACE)
target_include_directories(helmholtz_eos INTERFACE ${HELMHOLTZ_EOS_DIR}/include)
target_link_directories(helmholtz_eos INTERFACE ${HELMHOLTZ_EOS_DIR}/lib)
target_link_libraries(helmholtz_eos INTERFACE -lhelmholtzeos_c -lgfortran)
