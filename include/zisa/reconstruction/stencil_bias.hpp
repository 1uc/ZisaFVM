// SPDX-License-Identifier: MIT
// Copyright (c) 2021 ETH Zurich, Luc Grosheintz-Laval

#ifndef STENCIL_BIAS_H_D7PZ0
#define STENCIL_BIAS_H_D7PZ0

#include <string>

namespace zisa {

enum class StencilBias { central, one_sided };

StencilBias deduce_bias(const std::string &b);

} // namespace zisa

#endif /* end of include guard */
