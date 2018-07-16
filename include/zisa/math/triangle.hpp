/* A triangle consists of three vertices, edge-lengths and an area.
 */

#ifndef TRIANGLE_H_ONON1
#define TRIANGLE_H_ONON1

#include <zisa/math/cartesian.hpp>

namespace zisa {

class Triangle {
public:
  XY A;
  XY B;
  XY C;

  XY center;

  double a;
  double b;
  double c;

  double volume;

public:
  Triangle(const XY &A, const XY &B, const XY &C);
};

XY barycenter(const Triangle &tri);
double circum_radius(const Triangle &tri);
double avg_moment(const Triangle &tri, int x_deg, int y_deg, int quad_deg);
Triangle reference_triangle();

} // namespace zisa
#endif /* end of include guard */