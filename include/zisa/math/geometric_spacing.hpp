// SPDX-License-Identifier: MIT
// Copyright (c) 2021 ETH Zurich, Luc Grosheintz-Laval

#ifndef ZISA_GEOMETRIC_SPACING_DI8ON_HPP
#define ZISA_GEOMETRIC_SPACING_DI8ON_HPP

#include <zisa/config.hpp>
#include <zisa/memory/array.hpp>

namespace zisa {

array<double, 1>
geometric_spacing(double x_end, double growth_factor, int_t n_points);

}

#endif // ZISA_GEOMETRIC_SPACING_HPP
