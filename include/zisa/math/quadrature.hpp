/* Quadrature rules.
 */

#ifndef QUADRATURE_H_LP57B
#define QUADRATURE_H_LP57B

#include <map>

#include <zisa/config.hpp>
#include <zisa/grid/grid.hpp>
#include <zisa/math/cartesian.hpp>
#include <zisa/math/edge.hpp>
#include <zisa/math/edge_rule.hpp>
#include <zisa/math/triangle.hpp>
#include <zisa/math/triangular_rule.hpp>
#include <zisa/memory/array.hpp>

namespace zisa {

template <class QR, class F, class Domain>
auto quadrature(const QR &qr, const F &f, const Domain &domain)
    -> decltype(f(std::declval<XY>())) {

  using return_t = decltype(f(std::declval<XY>()));

  const auto &w = qr.weights;
  const auto &x = qr.points;

  // avoids assigning 'zero' to `ret`.
  return_t ret = static_cast<return_t>(w[0] * f(coord(domain, x[0])));

  // starts at --v
  for (int_t i = 1; i < qr.weights.size(); ++i) {
    ret = return_t(ret + w[i] * f(coord(domain, x[i])));
  }

  return return_t(volume(domain) * ret);
}

// -----------------
// -- Over Triangles

template <class F>
auto quadrature(const F &f, const Triangle &tri, int_t deg)
    -> decltype(f(std::declval<XY>())) {

  const auto &qr = cached_triangular_quadrature_rule(deg);
  return quadrature(qr, f, tri);
}

// -----------------
// -- Over Edges

template <class F>
auto quadrature(const F &f, const Edge &edge, int_t deg)
    -> decltype(f(std::declval<XY>())) {

  const auto &qr = cached_edge_quadrature_rule(deg);
  return quadrature(qr, f, edge);
}

// -----------------
// -- Over Grid

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

// -----------------
// -- average
template <class F, class Domain>
auto average(const F &f, const Domain &domain, int_t deg)
    -> decltype(f(std::declval<XY>())) {

  using return_t = decltype(f(std::declval<XY>()));

  auto vol = volume(domain);
  return return_t(1.0 / vol * quadrature(f, domain, deg));
}

} // namespace zisa

#endif /* end of include guard */
