#include <catch/catch.hpp>
#include <numeric>

#include <zisa/grid/grid.hpp>
#include <zisa/math/basic_functions.hpp>
#include <zisa/math/poly2d.hpp>
#include <zisa/math/quadrature.hpp>
#include <zisa/reconstruction/weno_ao.hpp>
#include <zisa/unit_test/math/basic_functions.hpp>

TEST_CASE("reconstruct smooth", "[weno_ao][math]") {
  auto f = [](const zisa::XY &x) {
    auto d = zisa::norm(x - zisa::XY{0.5, 0.5});
    double sigma = 0.2;
    return zisa::exp(-zisa::pow<2>(d / sigma));
  };

  auto qbar_local = zisa::array<double, 1>{zisa::shape_t<1>{20ul}};

  auto grid_names
      = std::vector<std::string>{"grids/convergence/unit_square_1.msh",
                                 "grids/convergence/unit_square_2.msh",
                                 "grids/convergence/unit_square_3.msh"};

  auto cases = std::vector<std::tuple<int, zisa::WENO_AO_Params>>{
      {1, {{{1}, {"c"}, {2.0}}, {1.0}}},
      {2, {{{2}, {"c"}, {2.0}}, {1.0}}},
      {3, {{{3}, {"c"}, {2.0}}, {1.0}}}};

  for (auto &[expected_rate, weno_ao_params] : cases) {
    std::vector<double> resolution;
    std::vector<double> l1_errors;
    std::vector<double> linf_errors;

    for (auto &&grid_name : grid_names) {
      auto grid = zisa::load_gmsh(grid_name);

      auto qbar = zisa::array<double, 1>(zisa::shape_t<1>{grid->n_cells});

      for (const auto &[i, tri] : triangles(*grid)) {
        qbar(i) = zisa::quadrature<4>(f, tri) / tri.volume;
      }

      double l1_err = 0.0;
      double linf_err = 0.0;
      double effective_volume = 0.0;

      for (const auto &[i, tri] : triangles(*grid)) {
        auto weno_ao = zisa::WENO_AO(grid, i, weno_ao_params);

        const auto &l2g = weno_ao.local2global();

        INFO(string_format("grid_name = '%s'", grid_name.c_str()));
        REQUIRE(l2g[0] == i);
        for (zisa::int_t i = 0; i < l2g.size(); ++i) {
          qbar_local(i) = qbar(l2g[i]);
        }

        auto p = weno_ao.reconstruct(qbar_local);

        auto diff = [tri = tri, &p, &f](const zisa::XY &x) {
          auto x_center = zisa::barycenter(tri);
          auto xx = zisa::XY((x - x_center) / zisa::circum_radius(tri));
          return zisa::abs(p(xx) - f(x));
        };

        auto err = zisa::quadrature<3>(diff, tri);

        l1_err += err;
        linf_err = zisa::max(linf_err, err);
        effective_volume += tri.volume;
      }

      resolution.push_back(zisa::largest_circum_radius(*grid));
      l1_errors.push_back(l1_err / effective_volume);
      linf_errors.push_back(linf_err);
    }

    auto rates = convergence_rates(resolution, l1_errors);

    for (zisa::int_t i = 0; i < rates.size(); ++i) {
      INFO(string_format(
          "[%d] (%e, %e) @ %e \n", i, l1_errors[i], rates[i], resolution[i]));
      CHECK(zisa::almost_equal(rates[i], expected_rate, 0.22));
    }
  }
}
