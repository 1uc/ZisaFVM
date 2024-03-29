// SPDX-License-Identifier: MIT
// Copyright (c) 2021 ETH Zurich, Luc Grosheintz-Laval

#include <zisa/testing/testing_framework.hpp>

#include <zisa/grid/grid.hpp>
#include <zisa/math/basic_functions.hpp>
#include <zisa/math/comparison.hpp>
#include <zisa/math/denormalized_rule.hpp>
#include <zisa/math/mathematical_constants.hpp>
#include <zisa/math/quadrature.hpp>

#include <zisa/unit_test/grid/test_grid_factory.hpp>
#include <zisa/unit_test/math/convergence_rates.hpp>

static void check_convergence(const std::vector<std::string> &grid_names,
                              double expected,
                              double atol,
                              zisa::int_t deg) {
  auto f = [](const zisa::XYZ &x) {
    // clang-format off
    return zisa::sin(0.5 * zisa::pi * x[0])
         + zisa::sin(0.5 * zisa::pi * x[1])
         + zisa::sin(0.5 * zisa::pi * x[2]);
    // clang-format on
  };

  std::vector<double> error;
  std::vector<double> resolution;

  for (auto &&grid_name : grid_names) {
    auto grid = zisa::load_grid(grid_name, deg);

    auto exact = grid->n_dims() * 2.0 / zisa::pi;
    auto approx = zisa::quadrature(f, *grid);
    error.push_back(zisa::abs(approx - exact));
    resolution.push_back(zisa::largest_circum_radius(*grid));
  }

  auto rates = convergence_rates(resolution, error);
  for (zisa::int_t i = 0; i < rates.size(); ++i) {
    if (error[i + 1] > atol) {
      INFO(string_format(
          "[%d] %e (%e), %e\n", i, rates[i], expected, error[i + 1]));
      REQUIRE(rates[i] > expected - 0.3);
    }
  }
}

static void
check_convergence_unit_square(double expected, double atol, zisa::int_t deg) {
  auto grid_names = zisa::TestGridFactory::unit_square();
  check_convergence(grid_names, expected, atol, deg);
}

static void
check_convergence_unit_cube(double expected, double atol, zisa::int_t deg) {
  auto grid_names = zisa::TestGridFactory::unit_cube();
  check_convergence(grid_names, expected, atol, deg);
}

TEST_CASE("Quadrature 2D; smooth f", "[math][quadrature][2d]") {
  SECTION("p == 1") { check_convergence_unit_square(2.0, 1e-14, 1); }
  SECTION("p == 2") { check_convergence_unit_square(3.0, 1e-14, 2); }
  SECTION("p == 3") { check_convergence_unit_square(4.0, 1e-14, 3); }
  SECTION("p == 4") { check_convergence_unit_square(5.0, 1e-14, 4); }
  SECTION("p == 5") { check_convergence_unit_square(6.0, 1e-14, 5); }
}

TEST_CASE("Quadrature 3D; smooth f", "[math][quadrature][3d]") {
  SECTION("p == 1") { check_convergence_unit_cube(2.0, 1e-14, 1); }
  SECTION("p == 2") { check_convergence_unit_cube(3.0, 1e-14, 2); }
  SECTION("p == 3") { check_convergence_unit_cube(4.0, 1e-14, 3); }
}

TEST_CASE("Quadrature; physical domain") {
  auto f = [](const zisa::XYZ &x) {
    return zisa::sin(0.5 * zisa::pi * x[0]) + zisa::sin(0.5 * zisa::pi * x[1]);
  };

  auto qr_ = zisa::cached_triangular_quadrature_rule(3);
  auto tri = zisa::reference_triangle();

  auto qr = zisa::denormalize(qr_, tri);

  auto fint_approx = zisa::quadrature(qr, f);
  auto fint_ref = zisa::quadrature(qr_, f, tri);
  REQUIRE(zisa::almost_equal(fint_approx, fint_ref, 1e-10));

  auto fbar_approx = zisa::average(qr, f);
  auto fbar_ref = zisa::average(qr_, f, tri);
  REQUIRE(zisa::almost_equal(fbar_approx, fbar_ref, 1e-10));
}