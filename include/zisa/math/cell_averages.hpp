#ifndef ZISA_CELL_AVERAGES_HPP
#define ZISA_CELL_AVERAGES_HPP

#include <zisa/math/quadrature.hpp>
#include <zisa/model/grid_variables.hpp>

namespace zisa {

template <class F>
GridVariables cell_averages(const F &f, const Grid &grid, int_t deg) {
  const auto &qr = cached_triangular_quadrature_rule(deg);

  int_t n_vars = decltype(f(std::declval<XYZ>()))::size();
  auto gvars = GridVariables(shape_t<2>{grid.n_cells, n_vars});
  for (const auto &[i, tri] : triangles(grid)) {
    gvars(i) = average(qr, f, tri);
  }

  return gvars;
}

}

#endif // ZISA_CELL_AVERAGES_HPP
