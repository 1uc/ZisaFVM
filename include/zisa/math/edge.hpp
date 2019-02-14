#ifndef EDGE_H_2L4RD
#define EDGE_H_2L4RD

#include <zisa/config.hpp>
#include <zisa/math/cartesian.hpp>

namespace zisa {

class Edge {
public:
  Edge(const XYZ &a, const XYZ &b)
      : a(a), b(b), n(rotate_right(normalize(b - a))), t(normalize(b - a)) {}

  Edge(const XYZ &a, const XYZ &b, const XYZ &n, const XYZ &t)
      : a(a), b(b), n(n), t(t) {}
  Edge(const Edge &other) = default;

  inline double volume() const { return zisa::norm(a - b); }

  inline const XYZ &normal() const { return n; }
  inline const XYZ &tangential() const { return t; }

  /// Quadrature reference space -> physical space.
  /** Input:
   *   x_rel : relative coordinate in [-1, 1].
   **/
  inline XYZ coord(double x_rel) const {
    // 0.5 * (a + b) + x_rel * (b - a)
    return XYZ((0.5 - 0.5 * x_rel) * a + (0.5 + 0.5 * x_rel) * b);
  }

private:
  XYZ a; // one end-point
  XYZ b; // other end-point

  XYZ n; // normal
  XYZ t; // tangential

private:
  friend XYZ unit_outward_normal(const Edge &edge, XYZ point_inside);
};

inline XYZ coord(const Edge &edge, double x_rel) { return edge.coord(x_rel); }
inline double volume(const Edge &edge) { return edge.volume(); }
inline XYZ unit_outward_normal(const Edge &edge, XYZ point_inside) {
  return XYZ(zisa::sign(zisa::dot(edge.normal(), edge.a - point_inside))
             * edge.normal());
}

} // namespace zisa

#endif
