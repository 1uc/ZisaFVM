#ifndef HYBRID_WENO_H_P7059
#define HYBRID_WENO_H_P7059

#include <string>
#include <vector>

#include <zisa/grid/grid.hpp>
#include <zisa/io/to_string.hpp>
#include <zisa/io/type_name.hpp>
#include <zisa/utils/indent_block.hpp>
#include <zisa/unit_test/math/basic_functions.hpp>
#include <zisa/unit_test/reconstruction/compute_convergence.hpp>

template <class RC>
void test_hybrid_weno_convergence(
    const std::vector<std::string> &grid_names,
    const std::tuple<double, double> &expected_rate,
    const zisa::HybridWENO_Params &params) {

  auto f = [](const zisa::XY &x) {
    auto d = zisa::norm(x - zisa::XY{0.5, 0.5});
    double sigma = 0.2;
    return zisa::exp(-zisa::pow<2>(d / sigma));
  };

  std::vector<double> resolution;
  std::vector<double> l1_errors;
  std::vector<double> linf_errors;

  for (auto &&grid_name : grid_names) {
    auto grid = zisa::load_gmsh(grid_name);

    auto rc = std::vector<RC>();
    for (const auto &[i, tri] : triangles(*grid)) {
      rc.push_back(RC(grid, i, params));
    }

    auto [l1_err, linf_err] = compute_errors(*grid, f, rc);

    resolution.push_back(zisa::largest_circum_radius(*grid));
    l1_errors.push_back(l1_err);
    linf_errors.push_back(linf_err);
  }

  auto rates = convergence_rates(resolution, l1_errors);

  for (zisa::int_t i = 0; i < rates.size(); ++i) {
    auto title = string_format("RC = %s", type_name<RC>().c_str());
    auto err_str = string_format("err[%d] = %e, rate[%d] = %e, res = %e",
                                 i,
                                 l1_errors[i],
                                 i,
                                 rates[i],
                                 resolution[i]);
    err_str = zisa::indent_block(1, err_str);

    auto desc_params = zisa::indent_block(1, zisa::to_string(params));
    INFO(string_format("%s\n%s\n%s", title.c_str(), err_str.c_str(), desc_params.c_str()));
    CHECK(zisa::is_inside_interval(rates[i], expected_rate));
  }
}

#endif /* end of include guard */
