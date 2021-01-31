cmake_minimum_required(VERSION 3.12 FATAL_ERROR)

project(Zisa VERSION 0.1)

include("cmake/has_build_type.cmake")
include("cmake/fast_debug.cmake")
include("cmake/native.cmake")
include("cmake/sanitizers.cmake")
include("cmake/gperftools.cmake")
include("cmake/clang_tidy.cmake")
include("cmake/warning_flags.cmake")



# -- Dependencies --------------------------------------------------------------
add_library(zisa_dependencies INTERFACE)

# experimental/filesystem
target_link_libraries(zisa_dependencies INTERFACE stdc++fs)

# Boost
find_package(Boost CONFIG REQUIRED program_options)
target_link_libraries(zisa_dependencies INTERFACE Boost::program_options)

# JSON
find_package(nlohmann_json MODULE REQUIRED)
target_link_libraries(zisa_dependencies INTERFACE nlohmann_json::nlohmann_json)

# Eigen
find_package(Eigen3 CONFIG REQUIRED)
target_link_libraries(zisa_dependencies INTERFACE Eigen3::Eigen)

# HDF5
find_package(HDF5 REQUIRED)
target_include_directories(zisa_dependencies INTERFACE ${HDF5_INCLUDE_DIRS})
target_link_libraries(zisa_dependencies INTERFACE ${HDF5_LIBRARIES})

# METIS
include("cmake/GetMetis.cmake")
target_link_libraries(zisa_dependencies INTERFACE metis)

# Helmholtz EOS
if(ZISA_HAS_HELMHOLTZ_EOS)
  include("cmake/GetHelmholtzEOS.cmake")
  target_compile_definitions(zisa_dependencies INTERFACE ZISA_HAS_HELMHOLTZ_EOS=1)
  target_link_libraries(zisa_dependencies INTERFACE helmholtz_eos)
endif()

# OpenGL
if(ZISA_HAS_OPENGL)
  find_package(GLEW REQUIRED)
  find_package(glfw3 REQUIRED)

  set(OpenGL_GL_PREFERENCE "GLVND")
  find_package(OpenGL REQUIRED)

  target_compile_definitions(zisa_dependencies INTERFACE ZISA_HAS_OPENGL=1)

  target_link_libraries(zisa_dependencies INTERFACE glfw GLEW::GLEW)
  target_link_libraries(zisa_dependencies INTERFACE ${OPENGL_LIBRARIES})
endif()

# Logging
if(ZISA_LOG_LEVEL)
  target_compile_definitions(zisa_dependencies INTERFACE ZISA_LOG_LEVEL=${ZISA_LOG_LEVEL})
endif()

# OpenMP
if(ZISA_HAS_OPENMP)
    list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
    find_package(OpenMP COMPONENTS CXX REQUIRED)

    target_compile_definitions(zisa_dependencies INTERFACE ZISA_HAS_OPENMP=1)
    target_link_libraries(zisa_dependencies INTERFACE OpenMP::OpenMP_CXX)
endif()

# MPI
if(ZISA_HAS_MPI)
    find_package(MPI REQUIRED)
    target_compile_definitions(zisa_dependencies INTERFACE ZISA_HAS_MPI=1)
    target_link_libraries(zisa_dependencies INTERFACE MPI::MPI_CXX)
endif()

# C++11 threads
find_package(Threads REQUIRED)
target_link_libraries(zisa_dependencies INTERFACE ${CMAKE_THREAD_LIBS_INIT})

# -- Zisa ----------------------------------------------------------------------
find_package(ZisaCore REQUIRED)
target_link_libraries(zisa_dependencies INTERFACE Zisa::core)

find_package(ZisaMemory REQUIRED)
target_link_libraries(zisa_dependencies INTERFACE Zisa::memory)

find_package(ZisaSFC REQUIRED)
target_link_libraries(zisa_dependencies INTERFACE Zisa::sfc)

find_package(ZisaTimeStepping REQUIRED)
target_link_libraries(zisa_dependencies INTERFACE Zisa::time_stepping)



# -- generic dependency --------------------------------------------------------
add_library(zisa_generic_obj OBJECT "")
set_target_properties(zisa_generic_obj PROPERTIES CXX_STANDARD 17)
target_include_directories(zisa_generic_obj PUBLIC include)
target_link_libraries(zisa_generic_obj PUBLIC zisa_dependencies)
target_link_libraries(zisa_generic_obj PRIVATE warning_flags)
target_link_libraries(zisa_generic_obj PUBLIC Zisa::core)

# -- MPI dependency ------------------------------------------------------------
if(ZISA_HAS_MPI)
  add_library(zisa_mpi_obj OBJECT "")
  set_target_properties(zisa_mpi_obj PROPERTIES CXX_STANDARD 17)
  target_include_directories(zisa_mpi_obj PUBLIC include)
  target_link_libraries(zisa_mpi_obj PUBLIC zisa_dependencies)
  target_link_libraries(zisa_mpi_obj PRIVATE warning_flags)
  target_link_libraries(zisa_mpi_obj PUBLIC Zisa::core)
endif()

# -- combined dependencies -----------------------------------------------------
add_library(zisa_obj INTERFACE)
target_sources(zisa_obj INTERFACE $<TARGET_OBJECTS:zisa_generic_obj>)

if(ZISA_HAS_MPI)
  target_sources(zisa_obj INTERFACE $<TARGET_OBJECTS:zisa_mpi_obj>)
endif()

target_include_directories(zisa_obj INTERFACE include)
target_link_libraries(zisa_obj INTERFACE zisa_dependencies)
target_link_libraries(zisa_obj INTERFACE Zisa::core)

# -- clang tidy ----------------------------------------------------------------
if(CLANG_TIDY_EXE)
  set_target_properties(zisa_obj PROPERTIES CXX_CLANG_TIDY "${DO_CLANG_TIDY}")
endif()

add_executable(zisa "")
set_target_properties(zisa PROPERTIES CXX_STANDARD 17)
target_link_libraries(zisa PUBLIC zisa_obj)

# -- Locate point --------------------------------------------------------------
add_executable(locate-point "")
set_target_properties(locate-point PROPERTIES CXX_STANDARD 17)
target_link_libraries(locate-point PUBLIC zisa_obj)

# -- OpenGL Demo ---------------------------------------------------------------
add_executable(opengl-demo "")
set_target_properties(opengl-demo PROPERTIES CXX_STANDARD 17)
target_link_libraries(opengl-demo PUBLIC zisa_obj)

# -- Something Demo ------------------------------------------------------------
add_executable(something-demo "")
set_target_properties(something-demo PROPERTIES CXX_STANDARD 17)
target_link_libraries(something-demo PUBLIC zisa_obj)

# -- Space Filling Curve Demo --------------------------------------------------
add_executable(space-filling-curve-demo "")
set_target_properties(space-filling-curve-demo PROPERTIES CXX_STANDARD 17)
target_link_libraries(space-filling-curve-demo PUBLIC zisa_obj)

# -- Domain Decomposition ------------------------------------------------------
add_executable(domain-decomposition "")
set_target_properties(domain-decomposition PROPERTIES CXX_STANDARD 17)
target_link_libraries(domain-decomposition PUBLIC zisa_obj metis)

# -- Renumber the Grid ---------------------------------------------------------
add_executable(renumber-grid "")
set_target_properties(renumber-grid PROPERTIES CXX_STANDARD 17)
target_link_libraries(renumber-grid PUBLIC zisa_obj metis)

# -- Test Reference solutions --------------------------------------------------
add_executable(test-distributed-reference "")
set_target_properties(test-distributed-reference PROPERTIES CXX_STANDARD 17)
target_link_libraries(test-distributed-reference PUBLIC zisa_obj)

# -- Tests ---------------------------------------------------------------------
find_package(Catch2 CONFIG REQUIRED)

add_executable(unit_tests "")
set_target_properties(unit_tests PROPERTIES CXX_STANDARD 17)
target_include_directories(unit_tests PUBLIC test)
target_link_libraries(unit_tests PUBLIC zisa_obj Catch2::Catch2)

enable_testing()
add_test(NAME UnitTests COMMAND unit_tests)

# -- Micro benchmarks ----------------------------------------------------------
if(ZISA_HAS_BENCHMARK)
  hunter_add_package(benchmark)
  find_package(benchmark CONFIG REQUIRED)

  add_executable(micro_benchmarks "")
  set_target_properties(micro_benchmarks PROPERTIES CXX_STANDARD 17)
  target_include_directories(micro_benchmarks PUBLIC benchmarks)
  target_link_libraries(micro_benchmarks PUBLIC zisa_obj benchmark::benchmark)
endif()

# -- Subdirectories ------------------------------------------------------------
add_subdirectory(src)
add_subdirectory(test)

if(ZISA_HAS_BENCHMARK)
  add_subdirectory(benchmarks)
endif()
