/* Quadrature rules.
 */

#ifndef QUADRATURE_H_LP57B
#define QUADRATURE_H_LP57B

#include <map>

#include <zisa/config.hpp>
#include <zisa/grid/grid.hpp>
#include <zisa/math/cartesian.hpp>
#include <zisa/math/cell.hpp>
#include <zisa/math/edge.hpp>
#include <zisa/math/edge_rule.hpp>
#include <zisa/math/face.hpp>
#include <zisa/math/triangle.hpp>
#include <zisa/math/triangular_rule.hpp>
#include <zisa/memory/array.hpp>

namespace zisa {

namespace detail {

struct UnitDomain {};
inline XYZ coord(const UnitDomain &, const XYZ &x) { return x; }
inline double volume(const UnitDomain &) { return 1.0; }

}

template <class QR, class F, class Domain>
auto quadrature(const QR &qr, const F &f, const Domain &domain)
    -> decltype(f(std::declval<XYZ>())) {

  using fx_t = decltype(f(std::declval<XYZ>()));

  const auto &w = qr.weights;
  const auto &x = qr.points;

  // avoids assigning 'zero' to `ret`.
  auto ret = fx_t(w[0] * f(coord(domain, x[0])));

  // starts at --v
  for (int_t i = 1; i < qr.weights.size(); ++i) {
    ret = fx_t(ret + w[i] * f(coord(domain, x[i])));
  }

  return fx_t(volume(domain) * ret);
}

template <class QR, class F>
auto quadrature(const QR &qr, const F &f) {
  return quadrature(qr, f, detail::UnitDomain{});
}

template <class QR, class F>
auto average(const QR &qr, const F &f) {
  using fx_t = decltype(f(std::declval<XYZ>()));
  return fx_t(quadrature(qr, f, detail::UnitDomain{}) / volume(qr));
}

// -----------------
// -- Over Triangle

template <class F>
auto quadrature(const F &f, const Triangle &tri, int_t deg)
    -> decltype(f(std::declval<XYZ>())) {

  const auto &qr = cached_triangular_quadrature_rule(deg);
  return quadrature(qr, f, tri);
}

// -----------------
// -- Over Edge

template <class F>
auto quadrature(const F &f, const Edge &edge, int_t deg)
    -> decltype(f(std::declval<XYZ>())) {

  const auto &qr = cached_edge_quadrature_rule(deg);
  return quadrature(qr, f, edge);
}

// -----------------
// -- Over Cell

template <class F>
auto quadrature(const Cell &cell, const F &f) {
  return quadrature(cell.qr, f);
}

template <class F>
auto average(const Cell &cell, const F &f) {
  return average(cell.qr, f);
}

// -----------------
// -- Over Face

template <class F>
auto quadrature(const Face &face, const F &f) {
  return quadrature(face.qr, f);
}

template <class F>
auto average(const Face &face, const F &f) {
  return average(face.qr, f);
}

// -----------------
// -- Over Grid

template <class F>
auto quadrature(const F &f, const Grid &grid)
    -> decltype(f(std::declval<XYZ>())) {

  // avoids assigning 'zero' to `ret`.
  auto ret = quadrature(grid.cells(0), f);

  // starts at --v
  for (int_t i : PlainIndexRange(1, grid.n_cells)) {
    ret = ret + quadrature(grid.cells(i), f);
  }

  return ret;
}

// -----------------
// -- average
template <class F, class Domain>
auto average(const F &f, const Domain &domain, int_t quad_deg)
    -> decltype(f(std::declval<XYZ>())) {

  using fx_t = decltype(f(std::declval<XYZ>()));
  return fx_t(quadrature(f, domain, quad_deg) / volume(domain));
}

template <class QR, class F, class Domain>
auto average(const QR &qr, const F &f, const Domain &domain)
    -> decltype(f(std::declval<XYZ>())) {

  using fx_t = decltype(f(std::declval<XYZ>()));
  return fx_t(quadrature(qr, f, domain) / volume(domain));
}

} // namespace zisa
#endif /* end of include guard */
