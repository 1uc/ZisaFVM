// SPDX-License-Identifier: MIT
// Copyright (c) 2021 ETH Zurich, Luc Grosheintz-Laval

#include <zisa/testing/testing_framework.hpp>

#include <random>

#include <zisa/io/colors.hpp>

TEST_CASE("Color Conversions") {
  zisa::RGBColor rgb{0.3f, 0.1f, 0.4f};
  auto xyz = zisa::rgb2xyz(rgb);

  float tol = 1e-4f;

  SECTION("RGB <-> XYZ") {
    auto approx = xyz2rgb(rgb2xyz(rgb));
    REQUIRE(almost_equal(rgb, approx, tol));
  }

  SECTION("XYZ <-> LAB") {
    auto approx = lab2xyz(xyz2lab(xyz));
    REQUIRE(almost_equal(xyz, approx, tol));
  }

  SECTION("RGB <-> LAB") {
    auto approx = lab2rgb(rgb2lab(rgb));
    REQUIRE(almost_equal(rgb, approx, tol));
  }
}
