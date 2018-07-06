#include <catch/catch.hpp>

#include <zisa/grid/grid.hpp>
#include <zisa/math/basic_functions.hpp>
#include <zisa/math/comparison.hpp>
#include <zisa/math/mathematical_constants.hpp>
#include <zisa/math/quadrature.hpp>

#include <zisa/unit_test/grid/test_grid_factory.hpp>
#include <zisa/unit_test/math/basic_functions.hpp>

template <int deg>
void check_convergence(double expected, double atol) {
  auto grid_names = zisa::TestGridFactory::unit_square();

  auto f = [](const zisa::XY &x) {
    return zisa::sin(0.5 * zisa::pi * x[0]) + zisa::sin(0.5 * zisa::pi * x[1]);
  };

  auto exact = 4.0 / zisa::pi;
  std::vector<double> error;
  std::vector<double> resolution;

  for (auto &&grid_name : grid_names) {
    auto grid = zisa::load_gmsh(grid_name);

    error.push_back(zisa::abs(zisa::quadrature<deg>(f, *grid) - exact));
    resolution.push_back(zisa::largest_circum_radius(*grid));
  }

  auto rates = convergence_rates(resolution, error);
  for (zisa::int_t i = 0; i < rates.size(); ++i) {
    if (error[i + 1] > atol) {
      INFO(string_format("[%d] %e (%e)\n", i, rates[i], expected));
      REQUIRE(std::abs(rates[i] - expected) < expected);
    }
  }
}

TEST_CASE("Quadrature; smooth f", "[math]") {
  SECTION("p == 1") { check_convergence<1>(2.0, 1e-14); }
  SECTION("p == 2") { check_convergence<2>(3.0, 1e-14); }
  SECTION("p == 3") { check_convergence<3>(4.0, 1e-14); }
  SECTION("p == 4") { check_convergence<4>(5.0, 1e-14); }
}
