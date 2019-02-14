/* Compute l1/linf errors. */

#ifndef COMPUTE_ERRORS_H_LJXFJ
#define COMPUTE_ERRORS_H_LJXFJ

#include <vector>

#include <zisa/grid/grid.hpp>
#include <zisa/math/quadrature.hpp>
#include <zisa/model/euler_variables.hpp>

namespace zisa {

template <class F, class RC>
std::tuple<double, double>
compute_errors(const Grid &grid, const F &f, const std::vector<RC> &rc) {

  constexpr int_t n_vars = 5;

  auto qbar_local = array<euler_var_t, 1>{shape_t<1>{1024ul}};
  auto qbar = array<euler_var_t, 1>(shape_t<1>{grid.n_cells});

  for (const auto &[i, tri] : triangles(grid)) {
    qbar(i) = zisa::average(f, tri, 4);
  }

  double l1_err = 0.0;
  double linf_err = 0.0;
  double volume = 0.0;

  for (const auto &[i, tri] : triangles(grid)) {
    const auto &l2g = rc[i].local2global();

    for (int_t ii = 0; ii < l2g.size(); ++ii) {
      assert(ii < qbar_local.size());
      for (int_t kk = 0; kk < n_vars; ++kk) {
        qbar_local(ii)[kk] = qbar(l2g[ii])[kk];
      }
    }

    auto p = rc[i].reconstruct(qbar_local);

    auto diff = [&p, &f](const XYZ &x) {
      return Cartesian<n_vars>{zisa::abs(p(x) - f(x))};
    };

    auto err = zisa::norm(quadrature(diff, tri, 3));

    l1_err += err;
    linf_err = zisa::max(linf_err, err);
  }

  return {l1_err, linf_err};
}
} // namespace zisa

#endif /* end of include guard */
