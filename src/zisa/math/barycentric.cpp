#include <zisa/math/barycentric.hpp>

namespace zisa {

Barycentric::Barycentric(double lambda1, double lambda2, double lambda3)
    : lambda{lambda1, lambda2, lambda3} {}

XYZ coord(const Triangle &tri, const Barycentric &lambda) {
  return XYZ(tri.A * lambda[0] + tri.B * lambda[1] + tri.C * lambda[2]);
}

std::ostream &operator<<(std::ostream &os, const Barycentric &bc) {
  os << string_format("[ %e, %e, %e]", bc[0], bc[1], bc[2]);
  return os;
}

} // namespace zisa
