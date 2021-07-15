// SPDX-License-Identifier: MIT
// Copyright (c) 2021 ETH Zurich, Luc Grosheintz-Laval

#include <zisa/reconstruction/stencil_bias.hpp>
#include <zisa/utils/logging.hpp>
#include <zisa/utils/string_format.hpp>

namespace zisa {

StencilBias deduce_bias(const std::string &b) {

  if (b == "c") {
    return StencilBias::central;
  } else if (b == "b") {
    return StencilBias::one_sided;
  }

  LOG_ERR(string_format("Missing case. [%s]", b.c_str()));
}

} // namespace zisa
