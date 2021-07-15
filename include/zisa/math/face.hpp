// SPDX-License-Identifier: MIT
// Copyright (c) 2021 ETH Zurich, Luc Grosheintz-Laval

#ifndef ZISA_FACE_HPP_JIDEW
#define ZISA_FACE_HPP_JIDEW

#include <zisa/config.hpp>

#include <utility>
#include <zisa/math/cartesian.hpp>
#include <zisa/math/denormalized_rule.hpp>

namespace zisa {
class Face {
public:
  DenormalizedRule qr;
  XYZ normal;
  std::pair<XYZ, XYZ> tangentials;

  Face() = default;
  Face(DenormalizedRule qr, const XYZ &n, const XYZ &t1, const XYZ &t2);
};

XYZ barycenter(const Face &face);

bool operator==(const Face &a, const Face &b);
bool operator!=(const Face &a, const Face &b);

double volume(const Face &face);

XYZ unit_outward_normal(const Face &face, const XYZ &point_inside);

}

#endif // ZISA_FACE_HPP
