// SPDX-License-Identifier: MIT
// Copyright (c) 2021 ETH Zurich, Luc Grosheintz-Laval

#ifndef ROLLING_CONVERGENCE_RATE_H_RBTXU
#define ROLLING_CONVERGENCE_RATE_H_RBTXU

#include <zisa/math/basic_functions.hpp>

namespace zisa {

template <class T>
class RollingConvergenceRate {
  static constexpr int max_values = 8;

public:
  RollingConvergenceRate(const T &atol) : atol(atol){};

  void push_back(const T &value) {
    values[i_end] = value;

    i_start = (n_values == max_values ? index(i_start + 1) : i_start);
    i_end = index(i_end + 1);
    n_values += 1;
  }

  bool is_converging(const T &expected_rate) const {
    if (zisa::all(values[index(i_end - 1)] <= atol)) {
      return true;
    }

    return zisa::all(rate(i_end - 2) >= expected_rate);
  }

private:
  auto rate(int i) const {
    const T &a = values[index(i - 1)];
    const T &b = values[index(i)];
    const T &c = values[index(i + 1)];

    return T(zisa::log(zisa::abs(c) / zisa::abs(b))
             / zisa::log(zisa::abs(b) / zisa::abs(a)));
  }

  int index(int i) const { return (i + max_values) % max_values; }

private:
  T atol;

  T values[max_values];
  int i_start = 0;
  int i_end = 0;
  int n_values = 0;
};

} // namespace zisa

#endif /* end of include guard */
