// SPDX-License-Identifier: MIT
// Copyright (c) 2021 ETH Zurich, Luc Grosheintz-Laval

#include <zisa/math/tetrahedral_rule.hpp>

#include <zisa/grid/grid.hpp>
#include <zisa/math/basic_functions.hpp>
#include <zisa/math/comparison.hpp>
#include <zisa/math/mathematical_constants.hpp>
#include <zisa/math/quadrature.hpp>
#include <zisa/testing/testing_framework.hpp>
#include <zisa/unit_test/grid/test_grid_factory.hpp>
#include <zisa/unit_test/math/convergence_rates.hpp>

TEST_CASE("TetrahedralRule; API", "[math]") {
  auto qr = zisa::make_tetrahedral_rule(3);
}

static void check_convergence(double expected, double atol, zisa::int_t deg) {
  auto grid_names = zisa::TestGridFactory::unit_cube();

  auto f = [](const zisa::XYZ &x) {
    // clang-format off
    return   zisa::sin(0.5  * zisa::pi * x[0])
           + zisa::cos(0.5  * zisa::pi * x[1])
           + zisa::sin(0.25 * zisa::pi * x[2]);
    // clang-format on
  };

  auto exact_A = 2.0 / zisa::pi;
  auto exact_B = 2.0 / zisa::pi;
  auto exact_C = 4.0 / zisa::pi * (zisa::cos(0.0) - zisa::cos(0.25 * zisa::pi));

  auto exact = exact_A + exact_B + exact_C;
  std::vector<double> error;
  std::vector<double> resolution;

  for (auto &&grid_name : grid_names) {
    auto grid = zisa::load_grid(grid_name, deg);

    error.push_back(zisa::abs(zisa::quadrature(f, *grid) - exact));
    resolution.push_back(zisa::largest_circum_radius(*grid));
  }

  auto rates = convergence_rates(resolution, error);
  for (zisa::int_t i = 0; i < rates.size(); ++i) {
    if (error[i + 1] > atol) {
      INFO(string_format(
          "[%d] %e != %e; %e\n", i, rates[i], expected, error[i + 1]));
      REQUIRE(rates[i] > expected - 0.3);
    }
  }
}

TEST_CASE("Quadrature tets; smooth f", "[math]") {
  SECTION("p == 1") { check_convergence(2.0, 1e-12, 1); }
  SECTION("p == 2") { check_convergence(3.0, 1e-12, 2); }
  SECTION("p == 3") { check_convergence(4.0, 2e-12, 3); }
}
