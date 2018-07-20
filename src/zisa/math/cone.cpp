#include <zisa/math/cone.hpp>
namespace zisa {

Cone::Cone(const XY &A, const XY &X, const XY &B)
    : x(X),
      dir(XY(normalize(normalize(A - X) + normalize(B - X)))),
      cos_angle(zisa::dot(normalize(A - X), dir)) {}

Cone::Cone(const XY &x, const XY &dir, double cos_angle)
    : x(x), dir(normalize(dir)), cos_angle(cos_angle) {}

bool Cone::is_inside(const XY &y) const {
  return zisa::dot(dir, normalize(y - x)) > cos_angle;
}

} // namespace zisa
