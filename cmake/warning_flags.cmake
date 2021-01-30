# -- Flags ---------------------------------------------------------------------
if(NOT TARGET warning_flags)
  add_library(warning_flags INTERFACE)

  if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    target_compile_options(warning_flags INTERFACE -Wall)
    target_compile_options(warning_flags INTERFACE -Wextra)
    target_compile_options(warning_flags INTERFACE -Wconversion)
  endif()

  if(ZISA_MAX_ERRORS)
    if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
      target_compile_options(warning_flags INTERFACE -fmax-errors=${ZISA_MAX_ERRORS})
    elseif(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
      target_compile_options(warning_flags INTERFACE -ferror-limit=${ZISA_MAX_ERRORS})
    endif()
  endif()
endif()


