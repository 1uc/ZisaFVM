// SPDX-License-Identifier: MIT
// Copyright (c) 2021 ETH Zurich, Luc Grosheintz-Laval

#include <zisa/math/face.hpp>

#include <zisa/math/edge.hpp>
#include <zisa/math/face_factory.hpp>
#include <zisa/math/triangle.hpp>

#include <zisa/testing/testing_framework.hpp>

TEST_CASE("Face; barycenter", "[math][face]") {

  auto quad_deg = zisa::int_t(3);
  double atol = 1e-10;

  SECTION("Edge") {
    auto v1 = zisa::XYZ{1.2, 3.0, 0.0};
    auto v2 = zisa::XYZ{4.2, 3.2, 4.0};

    auto edge = zisa::Edge(v1, v2);
    auto face = zisa::make_face(edge, quad_deg);

    auto approx = zisa::barycenter(face);
    auto exact = zisa::XYZ(0.5 * (v1 + v2));

    REQUIRE(zisa::almost_equal(approx, exact, atol));
  }

  SECTION("Edge") {
    auto v1 = zisa::XYZ{1.2, 3.0, 0.0};
    auto v2 = zisa::XYZ{4.2, 3.2, 4.0};
    auto v3 = zisa::XYZ{2.2, 3.2, 4.0};

    auto tri = zisa::Triangle(v1, v2, v3);
    auto face = zisa::make_face(tri, quad_deg);

    auto approx = zisa::barycenter(face);
    auto exact = zisa::XYZ((v1 + v2 + v3) / 3.0);

    REQUIRE(zisa::almost_equal(approx, exact, atol));
  }
}