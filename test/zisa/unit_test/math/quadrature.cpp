#include <zisa/testing/testing_framework.hpp>

#include <zisa/grid/grid.hpp>
#include <zisa/math/basic_functions.hpp>
#include <zisa/math/comparison.hpp>
#include <zisa/math/denormalized_rule.hpp>
#include <zisa/math/mathematical_constants.hpp>
#include <zisa/math/quadrature.hpp>

#include <zisa/unit_test/grid/test_grid_factory.hpp>
#include <zisa/unit_test/math/convergence_rates.hpp>

static void check_convergence(double expected, double atol, zisa::int_t deg) {
  auto grid_names = zisa::TestGridFactory::unit_square();

  auto f = [](const zisa::XYZ &x) {
    return zisa::sin(0.5 * zisa::pi * x[0]) + zisa::sin(0.5 * zisa::pi * x[1]);
  };

  auto exact = 4.0 / zisa::pi;
  std::vector<double> error;
  std::vector<double> resolution;

  for (auto &&grid_name : grid_names) {
    auto grid = zisa::load_gmsh(grid_name, deg);

    error.push_back(zisa::abs(zisa::quadrature(f, *grid) - exact));
    resolution.push_back(zisa::largest_circum_radius(*grid));
  }

  auto rates = convergence_rates(resolution, error);
  for (zisa::int_t i = 0; i < rates.size(); ++i) {
    if (error[i + 1] > atol) {
      INFO(string_format("[%d] %e (%e)\n", i, rates[i], expected));
      REQUIRE(rates[i] > expected - 0.3);
    }
  }
}

TEST_CASE("Quadrature; smooth f", "[math]") {
  SECTION("p == 1") { check_convergence(2.0, 1e-14, 1); }
  SECTION("p == 2") { check_convergence(3.0, 1e-14, 2); }
  SECTION("p == 3") { check_convergence(4.0, 1e-14, 3); }
  SECTION("p == 4") { check_convergence(5.0, 1e-14, 4); }
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