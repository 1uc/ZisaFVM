/* A triangle consists of three vertices, edge-lengths and an area.
 */

#ifndef TRIANGLE_H_ONON1
#define TRIANGLE_H_ONON1

#include <zisa/math/cartesian.hpp>
#include <zisa/math/edge.hpp>

namespace zisa {

class Triangle {
public:
  XYZ A;
  XYZ B;
  XYZ C;

  double a;
  double b;
  double c;

  double volume;

public:
  Triangle() = default;
  Triangle(const XYZ &A, const XYZ &B, const XYZ &C);
};

double volume(const Triangle &tri);
XYZ barycenter(const Triangle &tri);
double circum_radius(const Triangle &tri);
double characteristic_length(const Triangle &tri);
double inradius(const Triangle &tri);
double avg_moment(const Triangle &tri, int x_deg, int y_deg, int_t quad_deg);
bool is_inside(const Triangle &tri, const XYZ &x);
Triangle reference_triangle();

} // namespace zisa
#endif /* end of include guard */
