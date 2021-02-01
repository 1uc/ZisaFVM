# -- Flags ---------------------------------------------------------------------
add_library(fvm_warning_flags INTERFACE)
if(NOT TARGET warning_flags)
  if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    target_compile_options(fvm_warning_flags INTERFACE -Wall)
    target_compile_options(fvm_warning_flags INTERFACE -Wextra)
    target_compile_options(fvm_warning_flags INTERFACE -Wconversion)
  endif()

  if(ZISA_MAX_ERRORS)
    if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
      target_compile_options(fvm_warning_flags INTERFACE -fmax-errors=${ZISA_MAX_ERRORS})
    elseif(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
      target_compile_options(fvm_warning_flags INTERFACE -ferror-limit=${ZISA_MAX_ERRORS})
    endif()
  endif()
else()
  target_link_libraries(fvm_warning_flags INTERFACE warning_flags)
endif()


