/* Compute l1/linf errors. */

#ifndef COMPUTE_ERRORS_H_LJXFJ
#define COMPUTE_ERRORS_H_LJXFJ

#include <vector>
#include <zisa/grid/grid.hpp>
#include <zisa/math/quadrature.hpp>
#include <zisa/model/euler_variables.hpp>
#include <zisa/testing/testing_framework.hpp>

namespace zisa {

template <class F, class RC>
std::tuple<double, double>
compute_errors(const Grid &grid, const F &f, const array<RC, 1> &rc) {
  constexpr int_t n_vars = 5;

  auto qbar_local = array<euler_var_t, 1>{shape_t<1>{1024ul}};
  auto qbar = array<euler_var_t, 1>(shape_t<1>{grid.n_cells});
  auto polys = array<WENOPoly, 1>(shape_t<1>{16});
  auto rhs
      = array<double, 2, row_major>(shape_t<2>{1024ul, WENOPoly::n_vars()});

  for (const auto &[i, cell] : cells(grid)) {
    qbar(i) = zisa::average(cell.qr, f);
    //    CHECK(zisa::isreal(qbar(i)));
  }

  double l1_err = 0.0;
  double linf_err = 0.0;

  for (const auto &[i, cell] : cells(grid)) {
    const auto &l2g = rc[i].local2global();

    for (int_t ii = 0; ii < l2g.size(); ++ii) {
      assert(ii < qbar_local.size());

      for (int_t kk = 0; kk < n_vars; ++kk) {
        qbar_local(ii)[kk] = qbar(l2g[ii])[kk];
      }
    }

    auto p = rc[i].reconstruct(rhs, polys, qbar_local);

    auto diff = [&p, &f](const XYZ &x) {
      return Cartesian<n_vars>{zisa::abs(p(x) - f(x))};
    };

    if (!grid.cell_flags(i).ghost_cell) {
      auto err = zisa::norm(quadrature(cell.qr, diff));

      INFO(string_format("i = %d, err = %e, x = %s, ghost = %d",
                         i,
                         err,
                         format_as_list(grid.cell_centers[i]).c_str(),
                         int(grid.cell_flags[i].ghost_cell)));
      CHECK(zisa::isreal(err));

      l1_err += err;
      linf_err = zisa::max(linf_err, err);
    }
  }

  return {l1_err, linf_err};
}

} // namespace zisa
#endif /* end of include guard */
