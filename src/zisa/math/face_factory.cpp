// SPDX-License-Identifier: MIT
// Copyright (c) 2021 ETH Zurich, Luc Grosheintz-Laval

#include <zisa/math/face_factory.hpp>

#include <zisa/math/denormalized_rule.hpp>
#include <zisa/math/edge_rule.hpp>
#include <zisa/math/triangular_rule.hpp>

namespace zisa {
Face make_face(const Edge &edge, int_t quad_deg) {
  auto qr = denormalize(cached_edge_quadrature_rule(quad_deg), edge);

  const auto &[v0, v1] = edge.points;
  auto n = rotate_right(normalize(v1 - v0));
  auto t1 = XYZ(normalize(v1 - v0));
  auto t2 = XYZ(zisa::cross(n, t1));

  return Face(qr, n, t1, t2);
}

Face make_face(const Triangle &tri, int_t quad_deg) {
  auto qr = denormalize(cached_triangular_quadrature_rule(quad_deg), tri);

  const auto &v0 = tri.A;
  const auto &v1 = tri.B;
  const auto &v2 = tri.C;

  auto n = XYZ(normalize(zisa::cross(v1 - v0, v2 - v0)));
  auto t1 = XYZ(normalize(v1 - v0));
  auto t2 = XYZ(zisa::cross(n, t1));

  return Face(qr, n, t1, t2);
}
}