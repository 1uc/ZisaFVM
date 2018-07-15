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

double avg_moment(const Triangle &tri, int x_deg, int y_deg, int quad_deg) {
  auto center = barycenter(tri);

  auto f = [x_deg, y_deg, &center](const XY &xy) {
    auto x = XY(xy - center);
    return zisa::pow(x[0], x_deg) * zisa::pow(x[1], y_deg);
  };

  return 1.0 / tri.volume * quadrature(f, tri, quad_deg);
}

XY barycenter(const Triangle &tri) { return XY((tri.A + tri.B + tri.C) / 3.0); }

double circum_radius(const Triangle &tri) {
  auto c = barycenter(tri);
  return zisa::max(
      zisa::norm(tri.A - c), zisa::norm(tri.B - c), zisa::norm(tri.C - c));
}

} // namespace zisa
