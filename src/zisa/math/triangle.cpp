#include <zisa/math/basic_functions.hpp>
#include <zisa/math/quadrature.hpp>
#include <zisa/math/triangle.hpp>

namespace zisa {

Triangle::Triangle(const XY &A, const XY &B, const XY &C)
    : A(A),
      B(B),
      C(C),
      a(zisa::norm(B - C)),
      b(zisa::norm(A - C)),
      c(zisa::norm(A - B)),
      volume(herons_formula(a, b, c)) {}

double Triangle::avg_moment(int x_deg, int y_deg) {

  auto center = barycenter(*this);

  auto f = [this, x_deg, y_deg, &center](const Barycentric &bc) {
    auto x = XY(bc(*this) - center);
    return zisa::pow(x[0], x_deg) * zisa::pow(x[1], y_deg);
  };

  return 1.0 / volume * quadrature<3>(f, *this);
}

XY barycenter(const Triangle &tri) { return XY((tri.A + tri.B + tri.C) / 3.0); }

double circum_radius(const Triangle &tri) {
  auto c = barycenter(tri);
  return zisa::max(
      zisa::norm(tri.A - c), zisa::norm(tri.B - c), zisa::norm(tri.C - c));
}

} // namespace zisa
