#include <catch/catch.hpp>
#include <numeric>

#include <zisa/grid/grid.hpp>
#include <zisa/math/basic_functions.hpp>
#include <zisa/math/poly2d.hpp>
#include <zisa/math/quadrature.hpp>
#include <zisa/reconstruction/weno_ao.hpp>
#include <zisa/unit_test/math/basic_functions.hpp>

TEST_CASE("select stencil", "[weno_ao]") {

  auto grid = zisa::load_gmsh("grids/convergence/unit_square_2.msh");
  zisa::int_t i_cell = 42;
  zisa::int_t n_points = 12;

  SECTION("central_stencil") {
    const auto &x_cell = grid->cell_centers(i_cell);

    auto stencil = zisa::central_stencil(*grid, i_cell, n_points);
    REQUIRE(stencil.size() == n_points);

    auto dist = [&x_cell, &grid](zisa::int_t i) {
      auto x = grid->cell_centers(i);
      return zisa::norm(x - x_cell);
    };

    auto is_in_stencil = [&stencil](zisa::int_t i) {
      return std::find(stencil.begin(), stencil.end(), i) != stencil.end();
    };

    std::sort(
        stencil.begin(), stencil.end(), [&dist](zisa::int_t i, zisa::int_t j) {
          return dist(i) < dist(j);
        });

    auto r_max = dist(stencil.back());

    for (zisa::int_t i = 0; i < grid->cell_centers.shape(0); ++i) {
      REQUIRE(((dist(i) > r_max) || is_in_stencil(i)));
    }
  }

  SECTION("biased_stencil; k = 0") {
    zisa::int_t k = 0;
    auto stencil = zisa::biased_stencil(*grid, i_cell, k, n_points);
    REQUIRE(stencil.size() == n_points);
  }

  SECTION("biased_stencil; k = 1") {
    zisa::int_t k = 1;
    auto stencil = zisa::biased_stencil(*grid, i_cell, k, n_points);
    REQUIRE(stencil.size() == n_points);
  }

  SECTION("biased_stencil; k = 2") {
    zisa::int_t k = 2;
    auto stencil = zisa::biased_stencil(*grid, i_cell, k, n_points);
    REQUIRE(stencil.size() == n_points);
  }
}

TEST_CASE("deduce_max_order", "[weno_ao]") {

  double factor = 2.0;

  auto s = std::vector<zisa::int_t>(1);
  REQUIRE(zisa::deduce_max_order(s, factor) == 1);

  s.resize(2);
  REQUIRE(zisa::deduce_max_order(s, factor) == 1);

  s.resize(4);
  REQUIRE(zisa::deduce_max_order(s, factor) == 2);

  s.resize(9);
  REQUIRE(zisa::deduce_max_order(s, factor) == 2);

  s.resize(10);
  REQUIRE(zisa::deduce_max_order(s, factor) == 3);

  s.resize(17);
  REQUIRE(zisa::deduce_max_order(s, factor) == 3);

  s.resize(18);
  REQUIRE(zisa::deduce_max_order(s, factor) == 4);

  s.resize(27);
  REQUIRE(zisa::deduce_max_order(s, factor) == 4);

  s.resize(28);
  REQUIRE(zisa::deduce_max_order(s, factor) == 5);
}

TEST_CASE("reconstruct smooth", "[weno_ao][math]") {
  auto f = [](const zisa::XY &x) {
    return zisa::sin(zisa::pi * x[0] + 0.8 * x[1] * x[0] * x[0]);
  };

  auto qbar_local = zisa::array<double, 1>{zisa::shape_t<1>{20ul}};

  std::vector<double> resolution;
  std::vector<double> l1_errors;
  std::vector<double> linf_errors;

  auto grid_names
      = std::vector<std::string>{"grids/convergence/unit_square_0.msh",
                                 "grids/convergence/unit_square_1.msh",
                                 "grids/convergence/unit_square_2.msh",
                                 "grids/convergence/unit_square_3.msh"};

  for (auto &&grid_name : grid_names) {
    auto grid = zisa::load_gmsh(grid_name);

    auto qbar = zisa::array<double, 1>(zisa::shape_t<1>{grid->n_cells});

    for (auto [i, tri] : triangles(*grid)) {
      qbar(i) = zisa::quadrature<4>(f, tri);
    }

    double l1_err = 0.0;
    double linf_err = 0.0;

    for (auto [i, tri] : triangles(*grid)) {
      auto weno_ao = zisa::WENO_AO(grid, i);

      const auto &l2g = weno_ao.local2global();

      for (zisa::int_t i = 0; i < l2g.size(); ++i) {
        qbar_local(i) = qbar(l2g[i]);
      }

      auto p = weno_ao.reconstruct(qbar_local);

      auto f = [tri=tri, &p](const zisa::XY &x) {
        auto xx
            = zisa::XY((x - zisa::barycenter(tri)) / zisa::circum_radius(tri));
        return p(xx);
      };

      auto err = zisa::abs(zisa::quadrature<3>(f, tri) - qbar(i));

      l1_err += err * tri.volume;
      linf_err = zisa::max(linf_err, err);
    }

    resolution.push_back(zisa::largest_circum_radius(*grid));
    l1_errors.push_back(l1_err);
    linf_errors.push_back(linf_err);
  }

  auto rates = convergence_rates(resolution, l1_errors);

  for (zisa::int_t i = 0; i < rates.size(); ++i) {
    INFO(string_format("[%d] %e @ dx = %e \n", i, l1_errors[i], resolution[i]));
    CHECK(zisa::almost_equal(rates[i], 2.0, 0.3));
  }

  CHECK(l1_errors.back() < 4e-5);
}
