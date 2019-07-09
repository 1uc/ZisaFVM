#ifndef EDGE_H_2L4RD
#define EDGE_H_2L4RD

#include <zisa/config.hpp>
#include <zisa/math/cartesian.hpp>

namespace zisa {

class Edge {
public:
  XYZ points[2];

public:
  Edge(const XYZ &a, const XYZ &b) : points{a, b} {}
  Edge(const Edge &other) = default;

  inline double volume() const { return zisa::norm(points[1] - points[0]); }

  inline const XYZ &start_point() const { return points[0]; }
  inline const XYZ &end_point() const { return points[1]; }
};

inline XYZ coord(const Edge &edge, double x_rel) {
  const auto &[a, b] = edge.points;
  // 0.5 * (a + b) + x_rel * (b - a)
  return XYZ((0.5 - 0.5 * x_rel) * a + (0.5 + 0.5 * x_rel) * b);
}

inline double volume(const Edge &edge) { return edge.volume(); }

std::ostream &operator<<(std::ostream &os, const Edge &edge);

bool is_intersecting(const Edge &a, const Edge &b);

} // namespace zisa

#endif
