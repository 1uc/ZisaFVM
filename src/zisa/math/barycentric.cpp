// SPDX-License-Identifier: MIT
// Copyright (c) 2021 ETH Zurich, Luc Grosheintz-Laval

#include <zisa/math/barycentric.hpp>

namespace zisa {

Barycentric2D::Barycentric2D(double lambda1, double lambda2, double lambda3)
    : lambda{lambda1, lambda2, lambda3} {}

Barycentric2D::Barycentric2D(const Triangle &tri, const XYZ &x) {
  auto BA = XYZ(tri.B - tri.A);
  auto CA = XYZ(tri.C - tri.A);
  auto xA = XYZ(x - tri.A);

  double inv = 1.0 / (BA[0] * CA[1] - CA[0] * BA[1]);

  lambda[1] = (xA[0] * CA[1] - CA[0] * xA[1]) * inv;
  lambda[2] = (BA[0] * xA[1] - xA[0] * BA[1]) * inv;
  lambda[0] = 1.0 - lambda[1] - lambda[2];
}

XYZ coord(const Triangle &tri, const Barycentric2D &lambda) {
  return XYZ(tri.A * lambda[0] + tri.B * lambda[1] + tri.C * lambda[2]);
}

bool is_inside(const Barycentric2D &lambda) {
  return lambda[0] >= 0.0 && lambda[1] >= 0.0 && lambda[2] >= 0.0;
}

std::ostream &operator<<(std::ostream &os, const Barycentric2D &bc) {
  os << string_format("[ %e, %e, %e]", bc[0], bc[1], bc[2]);
  return os;
}

Barycentric3D::Barycentric3D(double lambda1,
                             double lambda2,
                             double lambda3,
                             double lambda4)
    : lambda{lambda1, lambda2, lambda3, lambda4} {}

Barycentric3D::Barycentric3D(const Tetrahedron &tet, const XYZ &x) {
  const auto &[v0, v1, v2, v3] = tet.points;

  auto dx = XYZ(x - v3);
  auto dv0 = XYZ(v0 - v3);
  auto dv1 = XYZ(v1 - v3);
  auto dv2 = XYZ(v2 - v3);

  auto [l0, l1, l2] = solve(dv0, dv1, dv2, dx);
  lambda[0] = l0;
  lambda[1] = l1;
  lambda[2] = l2;
  lambda[3] = 1.0 - (l0 + l1 + l2);
}

XYZ coord(const Tetrahedron &tet, const Barycentric3D &lambda) {
  return XYZ(tet.points[0] * lambda[0] + tet.points[1] * lambda[1]
             + tet.points[2] * lambda[2] + tet.points[3] * lambda[3]);
}

bool is_inside(const Barycentric3D &lambda) {
  return lambda[0] >= 0.0 && lambda[1] >= 0.0 && lambda[2] >= 0.0
         && lambda[3] >= 0.0;
}

std::ostream &operator<<(std::ostream &os, const Barycentric3D &bc) {
  os << string_format("[ %e, %e, %e]", bc[0], bc[1], bc[2]);
  return os;
}

} // namespace zisa
