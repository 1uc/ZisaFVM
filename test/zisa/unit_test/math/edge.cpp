#include <catch/catch.hpp>

#include <zisa/math/cartesian.hpp>
#include <zisa/math/edge.hpp>

using zisa::almost_equal;

TEST_CASE("Edge") {

  zisa::XY a{1.0, 2.0};
  zisa::XY b{2.0, 3.0};

  auto edge = zisa::Edge(a, b);
  double tol = 1e-12;

  REQUIRE(almost_equal(zisa::coord(edge, 0.0), zisa::XY(0.5 * (a + b)), tol));
  REQUIRE(almost_equal(zisa::coord(edge, -1.0), a, tol));
  REQUIRE(almost_equal(zisa::coord(edge, 1.0), b, tol));

  // |    b
  // |  a
  // |
  // --------
  REQUIRE(zisa::dot(edge.normal(), zisa::XY{1.0, -1.0}) > 0.0);
  REQUIRE(zisa::dot(edge.tangential(), b - a) > 0.0);
}
