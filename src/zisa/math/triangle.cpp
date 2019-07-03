#include <zisa/math/barycentric.hpp>
#include <zisa/math/basic_functions.hpp>
#include <zisa/math/quadrature.hpp>
#include <zisa/math/triangle.hpp>

namespace zisa {

Triangle::Triangle(const XYZ &A, const XYZ &B, const XYZ &C)
    : A(A),
      B(B),
      C(C),
      a(zisa::norm(B - C)),
      b(zisa::norm(A - C)),
      c(zisa::norm(A - B)),
      volume(herons_formula(a, b, c)) {}

double avg_moment(const Triangle &tri, int x_deg, int y_deg, int_t quad_deg) {
  auto center = barycenter(tri);

  auto f = [x_deg, y_deg, &center](const XYZ &x) {
    auto x_rel = XYZ(x - center);
    return zisa::pow(x_rel[0], x_deg) * zisa::pow(x_rel[1], y_deg);
  };

  return quadrature(f, tri, quad_deg) / tri.volume;
}

bool is_inside(const Triangle &tri, const XYZ &x) {
  return is_inside(Barycentric2D(tri, x));
}

double volume(const Triangle &tri) { return tri.volume; }

XYZ barycenter(const Triangle &tri) {
  return XYZ((tri.A + tri.B + tri.C) / 3.0);
}

double circum_radius(const Triangle &tri) {
  auto c = barycenter(tri);
  return zisa::max(
      zisa::norm(tri.A - c), zisa::norm(tri.B - c), zisa::norm(tri.C - c));
}

double inradius(const Triangle &tri) {
  double a = tri.a, b = tri.b, c = tri.c;
  return 0.5
         * zisa::sqrt((b + c - a) * (c + a - b) * (a + b - c) / (a + b + c));
}

double characteristic_length(const Triangle &tri) { return circum_radius(tri); }

Triangle reference_triangle() {
  return {{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}};
}

std::ostream &operator<<(std::ostream &os, const Triangle &tri) {
  os << "Triangle{" << tri.A << ", " << tri.B << ", " << tri.C << "}";
  return os;
}

} // namespace zisa
