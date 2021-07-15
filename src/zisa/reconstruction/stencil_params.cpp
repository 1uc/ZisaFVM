// SPDX-License-Identifier: MIT
// Copyright (c) 2021 ETH Zurich, Luc Grosheintz-Laval

#include <zisa/reconstruction/stencil_params.hpp>

namespace zisa {

StencilParams::StencilParams(int order, std::string bias, double overfit_factor)
    : order(order), bias(std::move(bias)), overfit_factor(overfit_factor) {}

} // namespace zisa
