#include <zisa/math/tetrahedron.hpp>

namespace zisa {
Tetrahedron::Tetrahedron(const XYZ &v0,
                         const XYZ &v1,
                         const XYZ &v2,
                         const XYZ &v3)
    : points{v0, v1, v2, v3} {}

double volume(const Tetrahedron &tetrahedron) {
  const auto &[v0, v1, v2, v3] = tetrahedron.points;

  auto [d11, d12, d13] = XYZ(v1 - v0);
  auto [d21, d22, d23] = XYZ(v2 - v0);
  auto [d31, d32, d33] = XYZ(v3 - v0);

  // clang-format off
  double vol = 1.0/6.0 * (d11*d22*d33 + d21*d32*d13 + d31*d12*d23
                        - d11*d32*d23 - d21*d12*d33 - d31*d22*d13);
  // clang-format on

  return zisa::abs(vol);
}

XYZ barycenter(const Tetrahedron &tetrahedron) {
  const auto &[v0, v1, v2, v3] = tetrahedron.points;
  return XYZ(0.25 * (v0 + v1 + v2 + v3));
}

}
