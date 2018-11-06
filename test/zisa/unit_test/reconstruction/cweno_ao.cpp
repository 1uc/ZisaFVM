#include <catch/catch.hpp>
#include <numeric>

#include <zisa/grid/grid.hpp>
#include <zisa/io/to_string.hpp>
#include <zisa/math/basic_functions.hpp>
#include <zisa/math/poly2d.hpp>
#include <zisa/math/quadrature.hpp>
#include <zisa/reconstruction/cweno_ao.hpp>
#include <zisa/unit_test/math/basic_functions.hpp>
#include <zisa/unit_test/reconstruction/compute_convergence.hpp>

TEST_CASE("CWENO_AO API", "[weno_ao][math][.]") {

  SECTION("compatibility with std::vector") {
    SECTION("push_back") {

      auto grid = zisa::load_gmsh("grids/small.msh");
      auto params = zisa::HybridWENO_Params({{{1}, {"c"}, {2.0}}, {1.0}});

      auto rc = std::vector<zisa::CWENO_AO>();
      for (const auto &[i, tri] : triangles(*grid)) {
        rc.push_back(zisa::CWENO_AO(grid, i, params));
      }

      REQUIRE(rc.size() == grid->n_cells);

      for (decltype(rc.size()) i = 0; i < rc.size(); ++i) {
        const auto &approx = rc[i];
        auto exact = zisa::CWENO_AO(grid, i, params);

        REQUIRE(approx == exact);
      }
    }
  }
}

// TEST_CASE("reconstruct smooth", "[weno_ao][math]") {
//   auto f = [](const zisa::XY &x) {
//     auto d = zisa::norm(x - zisa::XY{0.5, 0.5});
//     double sigma = 0.2;
//     return zisa::exp(-zisa::pow<2>(d / sigma));
//   };

//   auto grid_names
//       = std::vector<std::string>{"grids/convergence/unit_square_1.msh",
//                                  "grids/convergence/unit_square_2.msh",
//                                  "grids/convergence/unit_square_3.msh"};

//   using interval_t = std::tuple<double, double>;
//   auto cases = std::vector<std::tuple<interval_t, zisa::HybridWENO_Params>>{
//       {{0.8, 1.15}, {{{1}, {"c"}, {2.0}}, {1.0}}},
//       {{0.8, 1.15}, {{{1}, {"b"}, {2.0}}, {1.0}}},
//       {{1.8, 2.2}, {{{2}, {"c"}, {2.0}}, {1.0}}},
//       {{1.8, 2.2}, {{{2}, {"b"}, {2.0}}, {1.0}}},
//       {{2.8, 3.25}, {{{3}, {"c"}, {2.0}}, {1.0}}},
//       {{2.8, 3.25}, {{{3}, {"b"}, {2.0}}, {1.0}}},
//       {{3.8, 4.4}, {{{4}, {"c"}, {2.0}}, {1.0}}}};

//   // This expected not to converge high-order, because biases stencils at the
//   // boundary might not be large enough to be high-order. In such cases the
//   // order is reduced accordingly.
//   cases.push_back({{1.0, 4.1}, {{{4}, {"b"}, {2.0}}, {1.0}}});

//   cases.push_back({{2.8, 3.25},
//                    {{{2, 2, 2, 3}, {"b", "b", "b", "c"}, {1.5, 1.5, 1.5, 2.0}},
//                     {1.0, 1.0, 1.0, 100.0}}});

//   // The reason for the second order convergence is the linear weights. They
//   // allow too much pollution from the second order stencils.
//   cases.push_back({{1.8, 4.2},
//                    {{{4, 2, 2, 2}, {"c", "b", "b", "b"}, {2.0, 1.5, 1.5, 1.5}},
//                     {100.0, 1.0, 1.0, 1.0}}});

//   cases.push_back({{3.8, 4.4},
//                    {{{4, 2, 2, 2}, {"c", "b", "b", "b"}, {2.0, 1.5, 1.5, 1.5}},
//                     {1000.0, 1.0, 1.0, 1.0}}});

//   for (auto &[expected_rate, weno_ao_params] : cases) {
//     auto desc_params = zisa::to_string(weno_ao_params);
//     std::vector<double> resolution;
//     std::vector<double> l1_errors;
//     std::vector<double> linf_errors;

//     for (auto &&grid_name : grid_names) {
//       auto grid = zisa::load_gmsh(grid_name);

//       auto rc = std::vector<zisa::WENO_AO>();
//       for (const auto &[i, tri] : triangles(*grid)) {
//         rc.push_back(zisa::WENO_AO(grid, i, weno_ao_params));
//       }

//       auto [l1_err, linf_err] = compute_errors(*grid, f, rc);

//       resolution.push_back(zisa::largest_circum_radius(*grid));
//       l1_errors.push_back(l1_err);
//       linf_errors.push_back(linf_err);
//     }

//     auto rates = convergence_rates(resolution, l1_errors);

//     for (zisa::int_t i = 0; i < rates.size(); ++i) {
//       auto err_str = string_format("err[%d] = %e, rate[%d] = %e, res = %e",
//                                    i,
//                                    l1_errors[i],
//                                    i,
//                                    rates[i],
//                                    resolution[i]);

//       INFO(string_format("%s \n%s", err_str.c_str(), desc_params.c_str()));
//       CHECK(zisa::is_inside_interval(rates[i], expected_rate));
//     }
//   }
// }
