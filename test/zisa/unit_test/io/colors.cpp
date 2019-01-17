#include <zisa/testing/testing_framework.hpp>

#include <random>

#include <zisa/io/colors.hpp>

TEST_CASE("Color Conversions") {
  zisa::RGBColor rgb{0.3, 0.1, 0.4};
  auto xyz = zisa::rgb2xyz(rgb);

  float tol = 1e-4;

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
