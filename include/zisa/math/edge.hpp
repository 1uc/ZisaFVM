#ifndef EDGE_H_2L4RD
#define EDGE_H_2L4RD

#include <zisa/config.hpp>
#include <zisa/math/cartesian.hpp>

namespace zisa {

class Edge {
public:
  Edge(const XY &a, const XY &b)
      : a(a), b(b), n(rotate_right(normalize(b - a))), t(normalize(b - a)) {}

  Edge(const XY &a, const XY &b, const XY &n, const XY &t) : a(a), b(b), n(n), t(t) {}
  Edge(const Edge &other) = default;

  inline double volume() const { return zisa::norm(a - b); }

  inline const XY &normal() const { return n; }
  inline const XY &tangential() const { return t; }

  /// Quadrature reference space -> physical space.
  /** Input:
   *   x_rel : relative coordinate in [-1, 1].
   **/
  inline XY coord(double x_rel) const {
    // 0.5 * (a + b) + x_rel * (b - a)
    return XY((0.5 - 0.5 * x_rel) * a + (0.5 + 0.5 * x_rel) * b);
  }

private:
  XY a; // one end-point
  XY b; // other end-point

  XY n; // normal
  XY t; // tangential
};

inline XY coord(const Edge &edge, double x_rel) { return edge.coord(x_rel); }
inline double volume(const Edge &edge) { return edge.volume(); }

} // namespace zisa

#endif
