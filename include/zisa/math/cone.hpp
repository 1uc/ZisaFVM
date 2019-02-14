#ifndef CONE_H_CKYAW
#define CONE_H_CKYAW

#include <zisa/math/cartesian.hpp>

namespace zisa {

class Cone {
public:
  /// A 3D code with orgin X and A, B on opposing sides of the cone.
  Cone(const XYZ &A, const XYZ &X, const XYZ &B);
  Cone(const XYZ &x, const XYZ &dir, double cos_angle);
  bool is_inside(const XYZ &y) const;

private:
  XYZ x;
  XYZ dir;
  double cos_angle;
};

} // namespace zisa
#endif /* end of include guard */
