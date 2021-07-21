// SPDX-License-Identifier: MIT
// Copyright (c) 2021 ETH Zurich, Luc Grosheintz-Laval

#ifndef SCALAR_PLOT_H_HW2VU
#define SCALAR_PLOT_H_HW2VU

#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#include <zisa/config.hpp>
#include <zisa/io/color_map.hpp>
#include <zisa/io/colors.hpp>
#include <zisa/io/tri_plot.hpp>

namespace zisa {

class ScalarPlot {
public:
  ScalarPlot(const array<XYZ, 1> &vertices,
             const array<int_t, 2> &vertex_indices,
             const std::string &title);

  template <class F>
  void operator()(const F &f) {

    double min_ = std::numeric_limits<double>::max();
    double max_ = std::numeric_limits<double>::lowest();

    auto n_cells = values.size();
    for (int_t i = 0; i < n_cells; ++i) {
      values(i) = f(i);

      min_ = zisa::min(values(i), min_);
      max_ = zisa::max(values(i), max_);
    }

    std::tie(min_, max_) = safe_min_max(min_, max_);

    auto normalize = [min_, max_](double y) {
      return 2.0 * ((y - min_) / (max_ - min_)) - 1.0;
    };

    for (int_t i = 0; i < n_cells; ++i) {
      colors(3 * i) = color_map(normalize(values(i)));
      colors(3 * i + 1) = colors(3 * i);
      colors(3 * i + 2) = colors(3 * i);
    }

    tri_plot.draw(colors);
  }

private:
  std::pair<double, double> safe_min_max(double min_, double max_);

  TriPlot tri_plot;
  ColorMap color_map;

  array<double, 1> values;
  array<RGBColor, 1> colors;
};

} // namespace zisa
#endif
