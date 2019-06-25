#include <zisa/math/tetrahedron.hpp>

#include <zisa/testing/testing_framework.hpp>

TEST_CASE("Tetrahedron; volume", "[math]") {
  auto v0 = zisa::XYZ{0.0, 0.0, 0.0};
  auto v1 = zisa::XYZ{1.0, 0.0, 0.0};
  auto v2 = zisa::XYZ{0.0, 1.0, 0.0};
  auto v3 = zisa::XYZ{0.0, 0.0, 1.0};

  double atol = 1e-10;
  double vol_exact = 1.0 / 6.0;

  SECTION("v1") {
    auto tet = zisa::Tetrahedron{v0, v1, v2, v3};
    REQUIRE(zisa::almost_equal(volume(tet), vol_exact, atol));
  }

  SECTION("v2") {
    auto tet = zisa::Tetrahedron{v1, v0, v2, v3};
    REQUIRE(zisa::almost_equal(volume(tet), vol_exact, atol));
  }

  SECTION("v3") {
    auto tet = zisa::Tetrahedron{v1, v2, v3, v0};
    REQUIRE(zisa::almost_equal(volume(tet), vol_exact, atol));
  }

  SECTION("v3") {
    auto tet = zisa::Tetrahedron{v2, v1, v3, v0};
    REQUIRE(zisa::almost_equal(volume(tet), vol_exact, atol));
  }

  SECTION("v4") {
    auto tet = zisa::Tetrahedron{v0, zisa::XYZ{-v1}, v2, v3};
    REQUIRE(zisa::almost_equal(volume(tet), vol_exact, atol));
  }

  SECTION("v4") {
    auto tet = zisa::Tetrahedron{zisa::XYZ{-v0}, v1, zisa::XYZ{-v2}, v3};
    REQUIRE(zisa::almost_equal(volume(tet), vol_exact, atol));
  }
}