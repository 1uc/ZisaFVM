#include <zisa/testing/testing_framework.hpp>

#include <zisa/math/basic_functions.hpp>
#include <zisa/math/triangle.hpp>

TEST_CASE("is_inside", "[math]") {
  SECTION("reference triangle") {
    auto tri = zisa::reference_triangle();

    REQUIRE(is_inside(tri, zisa::XYZ{0.1, 0.1, 0.0}));
    REQUIRE(!is_inside(tri, zisa::XYZ{-0.1, 0.1, 0.0}));
    REQUIRE(!is_inside(tri, zisa::XYZ{1.1, 1.1, 0.0}));
  }

  SECTION("troublesome triangle") {
    auto tri = zisa::Triangle{{0.461024, 0.432542, 0.0},
                              {0.558699, 0.333827, 0.0},
                              {0.554918, 0.426801, 0.0}};
    auto x_center = zisa::XYZ{{0.52488, 0.397723, 0.0}};

    REQUIRE(is_inside(tri, x_center));
  }
}
