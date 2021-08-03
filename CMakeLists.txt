cmake_minimum_required(VERSION 3.14 FATAL_ERROR)

project(ZisaFVM VERSION 0.1)

if(ZISA_HAS_CUDA)
  enable_language(CUDA)
endif()

include("cmake/has_build_type.cmake")
include("cmake/fast_debug.cmake")
include("cmake/native.cmake")
include("cmake/sanitizers.cmake")
include("cmake/gperftools.cmake")
include("cmake/clang_tidy.cmake")
include("cmake/warning_flags.cmake")

# -- Dependencies --------------------------------------------------------------
add_library(fvm_dependencies INTERFACE)

# experimental/filesystem
target_link_libraries(fvm_dependencies INTERFACE stdc++fs)

# Boost
find_package(Boost REQUIRED program_options)
target_link_libraries(fvm_dependencies INTERFACE Boost::program_options)

# JSON
find_package(nlohmann_json REQUIRED)
target_link_libraries(fvm_dependencies INTERFACE nlohmann_json::nlohmann_json)

# Eigen
find_package(Eigen3 REQUIRED)
target_link_libraries(fvm_dependencies INTERFACE Eigen3::Eigen)

# HDF5
find_package(HDF5 REQUIRED)
target_include_directories(fvm_dependencies INTERFACE ${HDF5_INCLUDE_DIRS})
target_link_libraries(fvm_dependencies INTERFACE ${HDF5_LIBRARIES})

# METIS
if(ZISA_HAS_METIS)
  include("cmake/GetMetis.cmake")
  target_compile_definitions(fvm_dependencies INTERFACE ZISA_HAS_METIS=1)
  target_link_libraries(fvm_dependencies INTERFACE metis)
endif()

# Helmholtz EOS
if(ZISA_HAS_HELMHOLTZ_EOS)
  include("cmake/GetHelmholtzEOS.cmake")
  target_compile_definitions(fvm_dependencies INTERFACE ZISA_HAS_HELMHOLTZ_EOS=1)
  target_link_libraries(fvm_dependencies INTERFACE helmholtz_eos)
endif()

# OpenGL
if(ZISA_HAS_OPENGL)
  find_package(GLEW REQUIRED)
  find_package(glfw3 REQUIRED)

  set(OpenGL_GL_PREFERENCE "GLVND")
  find_package(OpenGL REQUIRED)

  target_compile_definitions(fvm_dependencies INTERFACE ZISA_HAS_OPENGL=1)

  target_link_libraries(fvm_dependencies INTERFACE glfw GLEW::GLEW)
  target_link_libraries(fvm_dependencies INTERFACE ${OPENGL_LIBRARIES})
endif()

# Logging
if(ZISA_LOG_LEVEL)
  target_compile_definitions(fvm_dependencies INTERFACE ZISA_LOG_LEVEL=${ZISA_LOG_LEVEL})
endif()

# OpenMP
if(ZISA_HAS_OPENMP)
    list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
    find_package(OpenMP COMPONENTS CXX REQUIRED)

    target_compile_definitions(fvm_dependencies INTERFACE ZISA_HAS_OPENMP=1)
    target_link_libraries(fvm_dependencies INTERFACE OpenMP::OpenMP_CXX)
endif()

# C++11 threads
find_package(Threads REQUIRED)
target_link_libraries(fvm_dependencies INTERFACE ${CMAKE_THREAD_LIBS_INIT})

# -- Zisa ----------------------------------------------------------------------
find_package(ZisaCore REQUIRED)
target_link_libraries(fvm_dependencies INTERFACE Zisa::core)

find_package(ZisaMemory REQUIRED)
target_link_libraries(fvm_dependencies INTERFACE Zisa::memory)

find_package(ZisaSFC REQUIRED)
target_link_libraries(fvm_dependencies INTERFACE Zisa::sfc)

find_package(ZisaTimeStepping REQUIRED)
target_link_libraries(fvm_dependencies INTERFACE Zisa::time_stepping)

if(ZISA_HAS_MPI)
  find_package(ZisaMPI REQUIRED)
  target_link_libraries(fvm_dependencies INTERFACE Zisa::mpi)
endif()

target_include_directories(fvm_dependencies INTERFACE
  $<BUILD_INTERFACE:${ZisaFVM_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>
)

# -- generic dependency --------------------------------------------------------
add_library(zisa_generic_obj OBJECT "")
set_target_properties(zisa_generic_obj PROPERTIES CXX_STANDARD 17)
set_target_properties(zisa_generic_obj PROPERTIES POSITION_INDEPENDENT_CODE ON)
target_link_libraries(zisa_generic_obj PUBLIC fvm_dependencies)
target_link_libraries(zisa_generic_obj PRIVATE fvm_warning_flags)

# -- MPI dependency ------------------------------------------------------------
if(ZISA_HAS_MPI)
  add_library(zisa_mpi_obj OBJECT "")
  set_target_properties(zisa_mpi_obj PROPERTIES CXX_STANDARD 17)
  set_target_properties(zisa_mpi_obj PROPERTIES POSITION_INDEPENDENT_CODE ON)
  target_link_libraries(zisa_mpi_obj PUBLIC fvm_dependencies)
  target_link_libraries(zisa_mpi_obj PRIVATE fvm_warning_flags)
endif()

# -- CUDA dependency ------------------------------------------------------------
if(ZISA_HAS_CUDA)
  add_library(zisa_cuda_obj OBJECT "")
  set_target_properties(zisa_cuda_obj PROPERTIES POSITION_INDEPENDENT_CODE ON)
  set_target_properties(zisa_cuda_obj PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
  set_target_properties(zisa_cuda_obj PROPERTIES CUDA_STANDARD 11)
  set_target_properties(zisa_cuda_obj PROPERTIES CXX_STANDARD 17)

  target_link_libraries(zisa_cuda_obj PUBLIC fvm_dependencies)
  target_link_libraries(zisa_cuda_obj PUBLIC Zisa::core)
  target_link_libraries(zisa_cuda_obj PUBLIC CUDA::cudart)
  target_link_libraries(zisa_cuda_obj PRIVATE fvm_warning_flags)
endif()

# -- combined dependencies -----------------------------------------------------
add_library(zisa_combined INTERFACE)
target_sources(zisa_combined INTERFACE $<TARGET_OBJECTS:zisa_generic_obj>)

if(ZISA_HAS_MPI)
  target_sources(zisa_combined INTERFACE $<TARGET_OBJECTS:zisa_mpi_obj>)
endif()

if(ZISA_HAS_CUDA)
  target_sources(zisa_combined INTERFACE $<TARGET_OBJECTS:zisa_cuda_obj>)
endif()


target_include_directories(zisa_combined INTERFACE include)
target_link_libraries(zisa_combined INTERFACE fvm_dependencies)

# # -- clang tidy ----------------------------------------------------------------
# if(CLANG_TIDY_EXE)
#   set_target_properties(zisa_generic PROPERTIES CXX_CLANG_TIDY "${DO_CLANG_TIDY}")
#   if(ZISA_HAS_MPI)
#   set_target_properties(zisa_mpi PROPERTIES CXX_CLANG_TIDY "${DO_CLANG_TIDY}")
#   endif()
# endif()

macro (zisa_generic_executable)
  message(${ARGV0})
  add_executable(${ARGV0} "")
  set_target_properties(${ARGV0} PROPERTIES CXX_STANDARD 17)
  set_target_properties(${ARGV0} PROPERTIES CUDA_STANDARD 11)
  set_target_properties(${ARGV0} PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
  target_link_libraries(${ARGV0} PUBLIC zisa_combined)
  target_link_libraries(${ARGV0} PRIVATE fvm_warning_flags)

  if(ZISA_HAS_CUDA)
    set_target_properties(${ARGV0} PROPERTIES CUDA_RESOLVE_DEVICE_SYMBOLS ON)
  endif()
endmacro()

zisa_generic_executable(zisa)

# -- Locate point --------------------------------------------------------------
zisa_generic_executable(locate-point)

# -- OpenGL Demo ---------------------------------------------------------------
zisa_generic_executable(opengl-demo)

# -- CUDA Demo ------------------------------------------------------------
zisa_generic_executable(cuda-demo)

# -- Something Demo ------------------------------------------------------------
zisa_generic_executable(something-demo)

# -- Space Filling Curve Demo --------------------------------------------------
zisa_generic_executable(space-filling-curve-demo)

# -- Domain Decomposition ------------------------------------------------------
zisa_generic_executable(domain-decomposition)

# -- Renumber the Grid ---------------------------------------------------------
zisa_generic_executable(renumber-grid)

# -- Test Reference solutions --------------------------------------------------
zisa_generic_executable(test-distributed-reference)

# -- Tests ---------------------------------------------------------------------
find_package(Catch2 REQUIRED)

zisa_generic_executable(unit_tests)
target_include_directories(unit_tests PUBLIC test)
target_link_libraries(unit_tests PUBLIC Catch2::Catch2)

enable_testing()
add_test(NAME UnitTests COMMAND unit_tests)

# -- Micro benchmarks ----------------------------------------------------------
if(ZISA_HAS_BENCHMARK)
  find_package(benchmark CONFIG REQUIRED)

  zisa_generic_executable(micro_benchmarks)
  target_include_directories(micro_benchmarks PUBLIC benchmarks)
  target_link_libraries(micro_benchmarks PUBLIC benchmark::benchmark)
endif()

# -- install -------------------------------------------------------------------
include(CMakePackageConfigHelpers)
write_basic_package_version_file(
  "${CMAKE_CURRENT_BINARY_DIR}/ZisaFVMConfigVersion.cmake"
  VERSION ${ZisaFVM_VERSION}
  COMPATIBILITY SameMajorVersion
  )

install(
  FILES
  "cmake/ZisaFVMConfig.cmake"
  "${CMAKE_CURRENT_BINARY_DIR}/ZisaFVMConfigVersion.cmake"
  DESTINATION lib/cmake/zisa
  )

install(
  TARGETS zisa
          locate-point
          renumber-grid
          domain-decomposition
  EXPORT ZisaFVMTargets
  )
include(GNUInstallDirs)
install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

install(EXPORT ZisaFVMTargets
  FILE ZisaFVMTargets.cmake
  NAMESPACE Zisa::
  DESTINATION lib/cmake/zisa
  )


# -- Subdirectories ------------------------------------------------------------
add_subdirectory(src)
add_subdirectory(test)

if(ZISA_HAS_BENCHMARK)
  add_subdirectory(benchmarks)
endif()
