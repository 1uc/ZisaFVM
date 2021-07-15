// SPDX-License-Identifier: MIT
// Copyright (c) 2021 ETH Zurich, Luc Grosheintz-Laval

#include <zisa/testing/testing_framework.hpp>

#include <zisa/math/quasi_newton.hpp>
#include <zisa/model/euler_variables.hpp>

TEST_CASE("quasi_newton; scalar", "[math][quasi_newton]") {
  double root = 1.2345;
  auto f = [root](double x) { return x * x - root * root; };
  auto df_inv
      = [](double x0) { return [x0](double fx) { return fx / (2.0 * x0); }; };
  double x0 = 1.9;
  double atol = 1e-10;

  auto [x_star, converged] = zisa::quasi_newton(f, df_inv, x0, atol);

  REQUIRE(converged);

  INFO(string_format("approx != exact | %e != %e", x_star, root));
  REQUIRE(zisa::almost_equal(x_star, root, 10.0 * atol));
}

TEST_CASE("quasi_newton; system", "[math][quasi_newton]") {
  auto theta = zisa::EnthalpyEntropy{1.0, 2.0};

  auto f = [](const zisa::EnthalpyEntropy &theta) {
    auto [h, K] = theta;

    return zisa::RhoE{h * K + h - K, K};
  };

  auto df_inv = [](const zisa::EnthalpyEntropy &theta0) {
    return [theta0](const zisa::RhoE &fx) {
      auto [h, K] = theta0;
      auto [rho, E] = fx;

      double inv_det = 1.0 / (K + 1.0);
      return zisa::EnthalpyEntropy{inv_det * (rho + (-h + 1.0) * E),
                                   inv_det * ((K + 1.0) * E)};
    };
  };

  auto x0 = zisa::EnthalpyEntropy{0.2, -0.4};
  auto atol = zisa::EnthalpyEntropy(1e-10);

  auto [x_star, converged] = zisa::quasi_newton(f, df_inv, x0, atol);

  REQUIRE(converged);

  auto root = zisa::EnthalpyEntropy(0.0);
  REQUIRE(zisa::almost_equal(x_star, root, 10.0 * atol[0]));
}
