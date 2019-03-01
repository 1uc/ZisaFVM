#include <zisa/math/barycentric.hpp>

namespace zisa {

Barycentric::Barycentric(double lambda1, double lambda2, double lambda3)
    : lambda{lambda1, lambda2, lambda3} {}

Barycentric::Barycentric(const Triangle &tri, const XYZ &x) {
  auto BA = XYZ(tri.B - tri.A);
  auto CA = XYZ(tri.C - tri.A);
  auto xA = XYZ(x - tri.A);

  double inv = 1.0 / (BA[0] * CA[1] - CA[0] * BA[1]);

  lambda[0] = (xA[0] * CA[1] - CA[0] * xA[1]) * inv;
  lambda[1] = (BA[0] * xA[1] - xA[0] * BA[1]) * inv;
  lambda[2] = 1.0 - lambda[0] - lambda[1];
}

XYZ coord(const Triangle &tri, const Barycentric &lambda) {
  return XYZ(tri.A * lambda[0] + tri.B * lambda[1] + tri.C * lambda[2]);
}

std::ostream &operator<<(std::ostream &os, const Barycentric &bc) {
  os << string_format("[ %e, %e, %e]", bc[0], bc[1], bc[2]);
  return os;
}

} // namespace zisa
