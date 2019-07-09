#include <zisa/testing/testing_framework.hpp>

#include <zisa/math/cartesian.hpp>
#include <zisa/math/edge.hpp>

using zisa::almost_equal;

TEST_CASE("Edge", "[math][edge]") {

  zisa::XYZ a{1.0, 2.0, 0.0};
  zisa::XYZ b{2.0, 3.0, 0.0};

  auto edge = zisa::Edge(a, b);
  double tol = 1e-12;

  REQUIRE(almost_equal(zisa::coord(edge, 0.0), zisa::XYZ(0.5 * (a + b)), tol));
  REQUIRE(almost_equal(zisa::coord(edge, -1.0), a, tol));
  REQUIRE(almost_equal(zisa::coord(edge, 1.0), b, tol));

  SECTION("is_intersection") {
    SECTION("edges") {
      auto edge_2
          = zisa::Edge{zisa::XYZ{3.5, 0.0, 0.0}, zisa::XYZ{0.0, 3.0, 0.0}};
      auto edge_3
          = zisa::Edge{zisa::XYZ{0.0, 0.0, 0.0}, zisa::XYZ{2.0, 1.0, 0.0}};

      REQUIRE(is_intersecting(edge, edge_2));
      REQUIRE(!is_intersecting(edge, edge_3));
    }

    SECTION("troublesome edges") {
      auto A = zisa::XYZ{0.461024, 0.432542, 0.0};
      auto B = zisa::XYZ{0.558699, 0.333827, 0.0};
      auto C = zisa::XYZ{0.554918, 0.426801, 0.0};

      auto e0 = zisa::Edge(A, B);
      auto e1 = zisa::Edge(B, C);
      auto e2 = zisa::Edge(C, A);

      auto x_inside = zisa::XYZ{0.52488, 0.397723, 0.0};
      auto x_outside = zisa::XYZ{0.6, 0.4, 0.0};

      auto crossing_edge = zisa::Edge(x_inside, x_outside);

      REQUIRE(!zisa::is_intersecting(e0, crossing_edge));
      REQUIRE(zisa::is_intersecting(e1, crossing_edge));
      REQUIRE(!zisa::is_intersecting(e2, crossing_edge));
    }
  }
}
