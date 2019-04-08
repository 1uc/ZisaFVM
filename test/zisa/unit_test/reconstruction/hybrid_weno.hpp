#ifndef HYBRID_WENO_H_P7059
#define HYBRID_WENO_H_P7059

#include <string>
#include <vector>

#include <zisa/grid/grid.hpp>
#include <zisa/reconstruction/global_reconstruction.hpp>
#include <zisa/unit_test/math/basic_functions.hpp>
#include <zisa/unit_test/reconstruction/compute_convergence.hpp>
#include <zisa/utils/indent_block.hpp>
#include <zisa/utils/to_string.hpp>
#include <zisa/utils/type_name.hpp>

namespace zisa {

template <class RC>
void test_hybrid_weno_convergence(
    const std::vector<std::string> &grid_names,
    const std::tuple<double, double> &expected_rate,
    const HybridWENOParams &params) {

  auto f = [](const XYZ &x) {
    auto d = zisa::norm(x - XYZ{0.5, 0.5, 0.0});
    double sigma = 0.2;
    double g = zisa::exp(-zisa::pow<2>(d / sigma));

    return euler_var_t{g, 0.5 * g, 0.25 * g, 0.4 * g, 0.35 * g};
  };

  std::vector<double> resolution;
  std::vector<double> l1_errors;
  std::vector<double> linf_errors;

  for (auto &&grid_name : grid_names) {
    auto grid = load_gmsh(grid_name);

    auto rc = std::vector<RC>();
    rc.reserve(grid->n_cells);

    for (const auto &[i, tri] : triangles(*grid)) {
      rc.push_back(RC(grid, i, params));
    }

    auto [l1_err, linf_err] = zisa::compute_errors(*grid, f, rc);

    resolution.push_back(largest_circum_radius(*grid));
    l1_errors.push_back(l1_err);
    linf_errors.push_back(linf_err);
  }

  auto rates = convergence_rates(resolution, l1_errors);

  for (int_t i = 0; i < rates.size(); ++i) {
    auto title = string_format("RC = %s", type_name<RC>().c_str());
    auto err_str = string_format("err[%d] = %e, rate[%d] = %e, res = %e",
                                 i,
                                 l1_errors[i],
                                 i,
                                 rates[i],
                                 resolution[i]);
    err_str = indent_block(1, err_str);

    auto desc_params = indent_block(1, zisa::to_string(params));
    INFO(string_format(
        "%s\n%s\n%s", title.c_str(), err_str.c_str(), desc_params.c_str()));
    CHECK(is_inside_interval(rates[i], expected_rate));
  }
}

template <class RC>
void test_hybrid_weno_stability(const std::vector<std::string> &grid_names,
                                const HybridWENOParams &params) {

  double tol = 5e-7;

  auto f = [](const XYZ &x) {
    auto d = zisa::norm(x - XYZ{0.5, 0.5, 0.0});
    return (d < 0.3 ? 10.0 : 1.0);
  };

  auto make_points = [](const Grid &grid, int_t i) {
    auto edge_rules = std::vector<EdgeRule>();
    auto volume_rules = std::vector<TriangularRule>();

    for (int_t deg = 1; deg <= 4; ++deg) {
      edge_rules.push_back(EdgeRule(deg));
      volume_rules.push_back(make_triangular_rule(deg));
    }

    auto points = std::vector<XYZ>();

    for (int_t k = 0; k < grid.max_neighbours; ++k) {
      auto edge = grid.edge(i, k);

      for (auto &&edge_rule : edge_rules) {
        for (auto &&x : edge_rule.points) {
          points.push_back(coord(edge, x));
        }
      }
    }

    auto tri = grid.triangle(i);
    for (auto &&volume_rule : volume_rules) {
      for (auto &&x : volume_rule.points) {
        points.push_back(coord(tri, x));
      }
    }

    return points;
  };

  for (auto &&grid_name : grid_names) {

    constexpr int_t n_vars = 5;

    auto grid = load_gmsh(grid_name);
    auto rc = EulerGlobalReconstruction<NoEquilibrium, RC>(
        grid, params, NoEquilibrium{});

    auto u = AllVariables({grid->n_cells, int_t(n_vars), int_t(0)});
    for (auto &&[i, tri] : triangles(*grid)) {
      for (int_t k = 0; k < n_vars; ++k) {
        u.cvars(i, k) = f(barycenter(tri));
      }
    }

    rc.compute(u);

    for (auto &&[i, tri] : triangles(*grid)) {
      auto points = make_points(*grid, i);
      for (auto &&x : points) {
        auto approx = rc(i)(x);
        auto exact = Cartesian<n_vars>(u.cvars(i));

        // clang-format off
        INFO(string_format("[%d] %s - %s = %s\n%s\n%s",
                           i,
                           zisa::to_string(approx).c_str(),
                           zisa::to_string(exact).c_str(),
                           zisa::to_string(zisa::Cartesian<n_vars>(approx - exact)).c_str(),
                           type_name<RC>().c_str(),
                           zisa::to_string(params).c_str()
                           ));
        // clang-format on

        REQUIRE(almost_equal(approx, exact, tol));
      }
    }
  }
}

} // namespace zisa
#endif /* end of include guard */
