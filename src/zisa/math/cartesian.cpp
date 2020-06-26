#include <zisa/math/cartesian.hpp>

namespace zisa {

Cartesian<3> solve(const Cartesian<3> &a,
                   const Cartesian<3> &b,
                   const Cartesian<3> &c,
                   const Cartesian<3> &x) {
  double idet = 1.0 / zisa::det(a, b, c);

  double A00 = b[1] * c[2] - b[2] * c[1];
  double A01 = a[1] * c[2] - a[2] * c[1];
  double A02 = a[1] * b[2] - a[2] * b[1];

  double A10 = b[0] * c[2] - b[2] * c[0];
  double A11 = a[0] * c[2] - a[2] * c[0];
  double A12 = a[0] * b[2] - a[2] * b[0];

  double A20 = b[0] * c[1] - b[1] * c[0];
  double A21 = a[0] * c[1] - a[1] * c[0];
  double A22 = a[0] * b[1] - a[1] * b[0];

  // clang-format off
  return Cartesian<3>({idet * ( A00 * x[0] - A10 * x[1] + A20 * x[2]),
                       idet * (-A01 * x[0] + A11 * x[1] - A21 * x[2]),
                       idet * ( A02 * x[0] - A12 * x[1] + A22 * x[2])});
  // clang-format on
}

}