// SPDX-License-Identifier: MIT
// Copyright (c) 2021 ETH Zurich, Luc Grosheintz-Laval

#include <zisa/math/cone.hpp>
#include <zisa/testing/testing_framework.hpp>

TEST_CASE("Cone", "[cone][math]") {

  SECTION("cone A") {
    auto c = zisa::TriangularCone(
        {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0});

    REQUIRE(c.is_inside({0.5, 0.5, 0.0}));
    REQUIRE(!c.is_inside({-0.5, 0.5, 0.0}));
    REQUIRE(!c.is_inside({0.5, -0.5, 0.0}));
    REQUIRE(!c.is_inside({-0.5, -0.5, 0.0}));
  }

  SECTION("cone B") {
    auto c = zisa::TriangularCone(
        {1.0, -2.0, 0.0}, {1.0, 1.0, 0.0}, {-0.0, -2.0, 0.0});

    REQUIRE(c.is_inside({0.5, 0.5, 0.0}));
    REQUIRE(c.is_inside({0.5, 0.5, 0.1}));
    REQUIRE(c.is_inside({-0.5, 0.5, 0.0}));
    REQUIRE(c.is_inside({-0.5, 0.5, 0.1}));

    REQUIRE(!c.is_inside({3.5, -0.5, 0.0}));
    REQUIRE(!c.is_inside({3.5, -0.5, 1.0}));
    REQUIRE(!c.is_inside({0.5, -4.5, 0.0}));
    REQUIRE(!c.is_inside({0.5, -4.5, 10.0}));
  }

  SECTION("cone C") {
    auto A = zisa::XYZ{0.0, 0.0, 0.0};
    auto B = zisa::XYZ{1.0, 0.0, 0.0};
    auto C = zisa::XYZ{1.0, 1.0, 0.0};
    auto D = zisa::XYZ{1.0, 0.0, 1.0};

    auto c = zisa::TetrahedralCone(A, B, C, D);

    REQUIRE(c.is_inside(zisa::XYZ(1.0 / 3.0 * (B + C + D))));
    REQUIRE(c.is_inside({3.0, 0.5, 0.1}));
    REQUIRE(!c.is_inside({-3.5, -0.5, 0.0}));
  }

  SECTION("cone D") {
    auto A = zisa::XYZ{0.4469, 0.7355, 0.0};
    auto B = zisa::XYZ{0.3894, 0.7428, 0.0};
    auto C = zisa::XYZ{0.4247, 0.6993, 0.0};

    auto c = zisa::TriangularCone(A, B, C);

    // points inside
    REQUIRE(c.is_inside({0.3727, 0.7158, 0.0}));
    REQUIRE(c.is_inside({0.3991, 0.7107, 0.0}));
    REQUIRE(c.is_inside({0.3630, 0.6848, 0.0}));

    // points outside
    REQUIRE(!c.is_inside({0.3582, 0.7625, 0.0}));
  }
}
