// SPDX-License-Identifier: MIT
// Copyright (c) 2021 ETH Zurich, Luc Grosheintz-Laval

#include <zisa/math/edge.hpp>

namespace zisa {

static bool ccw(const XYZ &A, const XYZ &B, const XYZ &C) {
  return (C[1] - A[1]) * (B[0] - A[0]) > (B[1] - A[1]) * (C[0] - A[0]);
}

bool is_intersecting(const zisa::Edge &a, const zisa::Edge &b) {

  const auto &A = a.start_point();
  const auto &B = a.end_point();
  const auto &C = b.start_point();
  const auto &D = b.end_point();

  return (ccw(A, C, D) != ccw(B, C, D)) and (ccw(A, B, C) != ccw(A, B, D));
}

std::ostream &operator<<(std::ostream &os, const Edge &edge) {
  os << "Edge(" << edge.start_point() << ", " << edge.end_point() << ")";
  return os;
}
}
