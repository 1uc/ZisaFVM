#include <zisa/math/cone.hpp>
namespace zisa {

Cone::Cone(const XYZ &A, const XYZ &X, const XYZ &B)
    : x(X),
      dir(XYZ(normalize(normalize(A - X) + normalize(B - X)))),
      cos_angle(zisa::dot(normalize(A - X), dir)) {}

Cone::Cone(const XYZ &x, const XYZ &dir, double cos_angle)
    : x(x), dir(normalize(dir)), cos_angle(cos_angle) {}

bool Cone::is_inside(const XYZ &y) const {
  return zisa::dot(dir, normalize(y - x)) > cos_angle;
}

} // namespace zisa
