#include <zisa/math/cone.hpp>
namespace zisa {

Cone::Cone(const XY &A, const XY &X, const XY &B)
    : x(X), dir(XY(normalize(A + B))), cos_angle(zisa::cos_angle(A, X, B)) {}

Cone::Cone(const XY &x, const XY &dir, double cos_angle)
    : x(x), dir(normalize(dir)), cos_angle(cos_angle) {}

bool Cone::is_inside(const XY &y) const {
  return zisa::dot(dir, y - x) > cos_angle;
}

} // namespace zisa
