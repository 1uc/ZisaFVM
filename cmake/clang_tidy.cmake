if(ZISA_HAS_CLANG_TIDY AND NOT DEFINED CLANG_TIDY_EXE)
  find_program(
          CLANG_TIDY_EXE
          NAMES "clang-tidy"
          DOC "Path to clang-tidy executable"
  )

  if(CLANG_TIDY_EXE)
    set(DO_CLANG_TIDY "${CLANG_TIDY_EXE}" "-checks=bugprone-*,cert-*,clang-analyser-*,mpi-*,performance-*,portability-*")
  else()
    message(FATAL_ERROR "Didn't find clang-tidy.")
  endif()
endif()

