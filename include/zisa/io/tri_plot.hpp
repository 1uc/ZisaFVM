// SPDX-License-Identifier: MIT
// Copyright (c) 2021 ETH Zurich, Luc Grosheintz-Laval

#ifndef TRI_PLOT_H_5752L
#define TRI_PLOT_H_5752L

#include <zisa/config.hpp>

#if ZISA_HAS_OPENGL != 0
#include <zisa/opengl/tri_plot.hpp>
#endif

#include <zisa/io/colors.hpp>
#include <zisa/math/cartesian.hpp>
#include <zisa/memory/array.hpp>

namespace zisa {

class TriPlot {
public:
  TriPlot(const array<XYZ, 1> &vertices,
          const array<int_t, 2> &vertex_indices,
          const std::string &title);

  void draw(const array<RGBColor, 1> &colors) const;

private:
#if ZISA_HAS_OPENGL != 0
  std::unique_ptr<opengl::TriPlot> tri_plot_;
#endif
};

} // namespace zisa

#endif
