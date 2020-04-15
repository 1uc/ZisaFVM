#include <zisa/testing/testing_framework.hpp>

#include <zisa/math/cartesian.hpp>
#include <zisa/model/euler_variables.hpp>

TEST_CASE("Cartesian; structured bindings", "[math]") {

  zisa::Cartesian<3> xyz{1.0, 2.0, 3.0};
  auto [x, y, z] = xyz;

  REQUIRE(x == xyz[0]);
  REQUIRE(y == xyz[1]);
  REQUIRE(z == xyz[2]);

  zisa::EnthalpyEntropy theta{1.0, 2.0};
  auto [h, K] = theta;

  REQUIRE(h == theta[0]);
  REQUIRE(K == theta[1]);

  zisa::euler_var_t u{1.0, 2.0, 3.0, 4.0, 5.0};
  auto [rho, mvx, mvy, mvz, E] = u;
}

TEST_CASE("Cartesian; comparison", "[math]") {
  zisa::Cartesian<3> x1{1.0, 2.0, 3.0};
  zisa::Cartesian<3> x2{1.0, 3.0, 3.0};

  REQUIRE(!zisa::all(x1 < x2));
  REQUIRE(zisa::all(x1 <= x2));
  REQUIRE(!zisa::all(x1 >= x2));
  REQUIRE(!zisa::all(x1 > x2));

  REQUIRE(zisa::any(x1 < x2));
  REQUIRE(!zisa::any(x1 > x2));
  REQUIRE(zisa::any(x1 >= x2));
  REQUIRE(zisa::any(x1 <= x2));
}
