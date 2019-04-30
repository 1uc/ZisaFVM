#ifndef ZISA_CELL_AVERAGES_HPP
#define ZISA_CELL_AVERAGES_HPP

#include <zisa/loops/for_each.hpp>
#include <zisa/math/quadrature.hpp>
#include <zisa/model/grid_variables.hpp>

namespace zisa {

template <class F>
GridVariables cell_averages(const F &f, const Grid &grid, int_t deg) {
  const auto &qr = cached_triangular_quadrature_rule(deg);

  int_t n_vars = decltype(f(std::declval<XYZ>()))::size();
  auto grid_vars = GridVariables(shape_t<2>{grid.n_cells, n_vars});

  zisa::for_each(triangles(grid),
                 [&grid_vars, &qr, &f](int_t i, const Triangle &tri) {
                   grid_vars(i) = average(qr, f, tri);
                 });

  return grid_vars;
}

}

#endif // ZISA_CELL_AVERAGES_HPP
