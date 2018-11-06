/* Compute l1/linf errors. */

#ifndef COMPUTE_ERRORS_H_LJXFJ
#define COMPUTE_ERRORS_H_LJXFJ

template <class F, class RC>
std::tuple<double, double>
compute_errors(const zisa::Grid &grid, const F &f, const std::vector<RC> &rc) {

  auto qbar_local = zisa::array<double, 1>{zisa::shape_t<1>{1024ul}};
  auto qbar = zisa::array<double, 1>(zisa::shape_t<1>{grid.n_cells});

  for (const auto &[i, tri] : triangles(grid)) {
    qbar(i) = zisa::quadrature<4>(f, tri) / tri.volume;
  }

  double l1_err = 0.0;
  double linf_err = 0.0;
  double volume = 0.0;

  for (const auto &[i, tri] : triangles(grid)) {
    const auto &l2g = rc[i].local2global();

    for (zisa::int_t ii = 0; ii < l2g.size(); ++ii) {
      assert(ii < qbar_local.size());
      qbar_local(ii) = qbar(l2g[ii]);
    }

    auto p = rc[i].reconstruct(qbar_local);

    auto diff = [tri = tri, &p, &f](const zisa::XY &x) {
      auto x_center = zisa::barycenter(tri);
      auto xx = zisa::XY((x - x_center) / zisa::circum_radius(tri));
      return zisa::abs(p(xx) - f(x));
    };

    auto err = zisa::quadrature<3>(diff, tri);

    l1_err += err;
    linf_err = zisa::max(linf_err, err);
  }

  return {l1_err, linf_err};
}

#endif /* end of include guard */