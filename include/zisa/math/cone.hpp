#ifndef CONE_H_CKYAW
#define CONE_H_CKYAW

#include<zisa/math/cartesian.hpp>

namespace zisa {

class Cone {
public:
  /// Cone AXB
  Cone(const XY &A, const XY &X, const XY &B);
  Cone(const XY &x, const XY &dir, double cos_angle);
  bool is_inside(const XY &y) const;

private:
  XY x;
  XY dir;
  double cos_angle;
};

} // namespace zisa
#endif /* end of include guard */
