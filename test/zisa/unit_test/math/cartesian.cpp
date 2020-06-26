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

TEST_CASE("Cartesian; solve 3x3", "[math][3x3]") {
  zisa::Cartesian<3> a0{1.0, 2.0, 3.0};
  zisa::Cartesian<3> a1{3.0, 0.0, 3.0};
  zisa::Cartesian<3> a2{0.0, -1.0, 0.0};
  zisa::Cartesian<3> x{2.32091, 0.2189, -0.2189};

  auto Ax = zisa::Cartesian<3>(x[0] * a0 + x[1] * a1 + x[2] * a2);
  auto x_approx = zisa::solve(a0, a1, a2, Ax);

  REQUIRE(zisa::almost_equal(x, x_approx, 1e-10));
}
