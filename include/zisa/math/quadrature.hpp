/* Quadrature rules.
 */

#ifndef QUADRATURE_H_LP57B
#define QUADRATURE_H_LP57B

#include <map>

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

  double operator[](int i) const { return lambda[i]; }

private:
  double lambda[3];
};

XY coord(const Triangle &tri, const Barycentric &x);

std::ostream &operator<<(std::ostream &os, const Barycentric &bc);

struct QuadratureRule {
  array<double, 1> weights;
  array<Barycentric, 1> points;

  QuadratureRule(int_t n_points);
  QuadratureRule(const QuadratureRule &qr) = default;
  QuadratureRule(QuadratureRule &&qr) = default;
};

QuadratureRule make_triangular_rule(int_t deg);

const QuadratureRule &cached_triangular_quadrature_rule(int_t deg);

template <class F>
auto quadrature(const QuadratureRule &qr, const F &f, const Triangle &tri)
    -> decltype(f(std::declval<XY>())) {

  const auto &w = qr.weights;
  const auto &x = qr.points;

  // avoids assigning 'zero' to `ret`.
  decltype(f(coord(tri, x[0]))) ret = w[0] * f(coord(tri, x[0]));

  // starts at --v
  for (int_t i = 1; i < qr.weights.size(); ++i) {
    ret = ret + w[i] * f(coord(tri, x[i]));
  }

  return tri.volume * ret;
}

template <class F>
auto quadrature(const F &f, const Triangle &tri, int_t deg)
    -> decltype(f(std::declval<XY>())) {

  const auto &qr = cached_triangular_quadrature_rule(deg);
  return quadrature(qr, f, tri);
}

template <class F>
auto quadrature(const F &f, const Grid &grid, int_t deg)
    -> decltype(f(std::declval<XY>())) {

  const auto &qr = cached_triangular_quadrature_rule(deg);

  auto tri = grid.triangle(0);

  // avoids assigning 'zero' to `ret`.
  auto ret = quadrature(f, tri, deg);

  // starts at --v
  for (int_t i = 1; i < grid.n_cells; ++i) {
    tri = grid.triangle(i);
    ret = ret + quadrature(qr, f, tri);
  }

  return ret;
}

} // namespace zisa

#endif /* end of include guard */
