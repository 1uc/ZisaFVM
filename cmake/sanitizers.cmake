has_build_type(ASAN)
if(NOT ZISA_HAS_DEFINED_ASAN)
  SET(CMAKE_CXX_FLAGS_ASAN "-O2 -g -fsanitize=address -fno-omit-frame-pointer")
endif()

has_build_type(UBSAN)
if(NOT ZISA_HAS_DEFINED_UBSAN)
  SET(CMAKE_CXX_FLAGS_UBSAN "-O2 -g -fsanitize=undefined")
  SET(CMAKE_EXE_LINKER_FLAGS_UBSAN "-fsanitize=undefined")
endif()
