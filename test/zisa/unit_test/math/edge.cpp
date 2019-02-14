#include <zisa/testing/testing_framework.hpp>

#include <zisa/math/cartesian.hpp>
#include <zisa/math/edge.hpp>

using zisa::almost_equal;

TEST_CASE("Edge") {

  zisa::XYZ a{1.0, 2.0, 0.0};
  zisa::XYZ b{2.0, 3.0, 0.0};

  auto edge = zisa::Edge(a, b);
  double tol = 1e-12;

  REQUIRE(almost_equal(zisa::coord(edge, 0.0), zisa::XYZ(0.5 * (a + b)), tol));
  REQUIRE(almost_equal(zisa::coord(edge, -1.0), a, tol));
  REQUIRE(almost_equal(zisa::coord(edge, 1.0), b, tol));

  // |    b
  // |  a
  // |
  // --------
  REQUIRE(zisa::dot(edge.normal(), zisa::XYZ{1.0, -1.0, 0.0}) > 0.0);
  REQUIRE(zisa::dot(edge.tangential(), b - a) > 0.0);
}
