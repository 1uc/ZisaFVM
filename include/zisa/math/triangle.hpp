// SPDX-License-Identifier: MIT
// Copyright (c) 2021 ETH Zurich, Luc Grosheintz-Laval

/* A triangle consists of three vertices, edge-lengths and an area.
 */

#ifndef TRIANGLE_H_ONON1
#define TRIANGLE_H_ONON1

#include <zisa/math/cartesian.hpp>
#include <zisa/math/edge.hpp>

namespace zisa {

class Triangle {
public:
  XYZ A;
  XYZ B;
  XYZ C;

  double a;
  double b;
  double c;

  double volume;

public:
  Triangle() = default;
  Triangle(const XYZ &A, const XYZ &B, const XYZ &C);
};

std::ostream &operator<<(std::ostream &os, const Triangle &tri);

double volume(const Triangle &tri);
XYZ barycenter(const Triangle &tri);
double circum_radius(const Triangle &tri);
double characteristic_length(const Triangle &tri);
double inradius(const Triangle &tri);
bool is_inside(const Triangle &tri, const XYZ &x);
Triangle reference_triangle();

} // namespace zisa
#endif /* end of include guard */
