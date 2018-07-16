/* Quadrature rules.
 */

#ifndef QUADRATURE_H_LP57B
#define QUADRATURE_H_LP57B

#include <zisa/config.hpp>
#include <zisa/grid/grid.hpp>
#include <zisa/math/cartesian.hpp>
#include <zisa/math/triangle.hpp>
#include <zisa/memory/array.hpp>

namespace zisa {

class Barycentric {
public:
  Barycentric() = default;
  Barycentric(double lambda1, double lambda2, double lambda3);

  XY operator()(const Triangle &tri) const;

  double operator[](int i) const { return lambda[i]; }

private:
  double lambda[3];
};

std::ostream &operator<<(std::ostream &os, const Barycentric &bc);

struct QuadratureRule {
  array<double, 1> weights;
  array<Barycentric, 1> points;

  QuadratureRule(int_t n_points);
};

QuadratureRule make_triangular_rule(int n_points);

template <int n_points>
const QuadratureRule &cached_triangular_quadrature_rule() {
  static auto qr = make_triangular_rule(n_points);
  return qr;
}

template <int deg, class F>
auto quadrature(const F &f, const Triangle &tri)
    -> decltype(f(std::declval<XY>())) {

  const auto &qr = cached_triangular_quadrature_rule<deg>();
  const auto &w = qr.weights;
  const auto &x = qr.points;

  decltype(f(x[0](tri))) ret = w[0] * f(x[0](tri));

  for (int_t i = 1; i < qr.weights.size(); ++i) {
    ret = ret + w[i] * f(x[i](tri));
  }

  return tri.volume * ret;
}

template <class F>
auto quadrature(const F &f, const Triangle &tri, int deg)
    -> decltype(f(std::declval<XY>())) {

  if (deg == 1) {
    return quadrature<1>(f, tri);
  } else if (deg == 2) {
    return quadrature<2>(f, tri);
  } else if (deg == 3) {
    return quadrature<3>(f, tri);
  } else if (deg == 4) {
    return quadrature<4>(f, tri);
  } else {
    LOG_ERR("Implement the remaining quadrature rules.");
  }
}

template <int deg, class F>
auto quadrature(const F &f, const Grid &grid)
    -> decltype(f(std::declval<XY>())) {

  auto tri = grid.triangle(0);

  auto ret = quadrature<deg>(f, tri);

  for (int_t i = 1; i < grid.n_cells; ++i) {
    tri = grid.triangle(i);
    ret = ret + quadrature<deg>(f, tri);
  }

  return ret;
}

} // namespace zisa

#endif /* end of include guard */
