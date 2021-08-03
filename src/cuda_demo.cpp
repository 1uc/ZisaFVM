// SPDX-License-Identifier: MIT
// Copyright (c) 2021 ETH Zurich, Luc Grosheintz-Laval

#include <zisa/config.hpp>
#include <zisa/cuda/hello_world.hpp>

int main() {
#if ZISA_HAS_CUDA > 0

  zisa::hello_world();

#else
  LOG_ERR("Compiled without CUDA support.");
#endif
}
