// SPDX-License-Identifier: MIT
// Copyright (c) 2021 ETH Zurich, Luc Grosheintz-Laval

#include <zisa/math/basic_functions.hpp>
#include <zisa/math/newton.hpp>

namespace zisa {

double newton(const std::function<double(double)> &f,
              const std::function<double(double)> &df,
              double x_guess,
              double atol,
              int max_iter) {

  return newton(
      [&f, &df](double x) {
        return std::pair<double, double>{f(x), df(x)};
      },
      x_guess,
      atol,
      max_iter);
}

}