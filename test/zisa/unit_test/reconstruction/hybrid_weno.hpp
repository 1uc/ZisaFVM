#ifndef HYBRID_WENO_H_P7059
#define HYBRID_WENO_H_P7059

#include <string>
#include <vector>

#include <zisa/grid/grid.hpp>
#include <zisa/io/hdf5_serial_writer.hpp>
#include <zisa/math/tetrahedral_rule.hpp>
#include <zisa/memory/array_stencil_family.hpp>
#include <zisa/reconstruction/global_reconstruction.hpp>
#include <zisa/unit_test/math/convergence_rates.hpp>
#include <zisa/unit_test/reconstruction/compute_convergence.hpp>
#include <zisa/utils/indent_block.hpp>
#include <zisa/utils/join.hpp>
#include <zisa/utils/to_string.hpp>
#include <zisa/utils/type_name.hpp>

namespace zisa {

inline void save_grid_with_stencil(const std::string &filename,
                                   const Grid &grid,
                                   const array<StencilFamily, 1> &stencils) {
  auto writer = HDF5SerialWriter(filename);
  save(writer, grid);
  save(writer, stencils, "stencils");
}

inline std::pair<std::shared_ptr<Grid>, array<StencilFamily, 1>>
load_grid_with_stencils(const std::vector<std::string> &grid_names,
                        int grid_level,
                        const HybridWENOParams &params,
                        int quad_deg) {

  auto grid_name = grid_names[grid_level];
  auto grid = load_grid(grid_name, quad_deg);

  mask_ghost_cells(*grid, [](const Grid &grid, int_t i) {
    auto x = grid.cell_centers(i);
    return zisa::any(x < Cartesian<3>(0.0)) || zisa::any(x > Cartesian<3>(1.0));
  });
  auto n_cells = grid->n_cells;

  auto sf = array<StencilFamily, 1>(shape_t<1>{n_cells});
  for (auto i : cell_indices(*grid)) {
    sf[i] = StencilFamily(*grid, i, params.stencil_family_params);
  }

  return {grid, std::move(sf)};
}

template <class RC>
void test_hybrid_weno_valid_stencil(
    const std::vector<std::string> &grid_names,
    const std::tuple<double, double> &expected_rate,
    const HybridWENOParams &params) {

  auto quad_deg = 4;

  for (int grid_level = 0; grid_level < grid_names.size(); ++grid_level) {
    auto [grid, sf]
        = load_grid_with_stencils(grid_names, grid_level, params, quad_deg);

    auto debug_filename = string_format(
        "__unit_tests-weno_ao-stencil-%d--%d.h5", grid_level, rand());
    save_grid_with_stencil(debug_filename, *grid, sf);

    for (auto i : cell_indices(*grid)) {
      if (!grid->cell_flags[i].ghost_cell) {
        const auto &orders = params.stencil_family_params.orders;
        for (int_t k = 0; k < orders.size(); ++k) {

          auto title = string_format("RC = %s", type_name<RC>().c_str());
          auto desc_params = indent_block(1, zisa::to_string(params));
          auto info = indent_block(
              1,
              string_format("(%d, %d) order = %d, size = %d, x = %s",
                            i,
                            k,
                            sf[i][k].order(),
                            sf[i][k].local().size(),
                            format_as_list(grid->cell_centers[i]).c_str()));
          INFO(string_format("%s\n%s\n%s\ndebug_file = %s",
                             title.c_str(),
                             desc_params.c_str(),
                             info.c_str(),
                             debug_filename.c_str()));
          REQUIRE(sf[i][k].order() == orders[k]);
        }
      }
    }
  }
}

template <class RC>
void test_hybrid_weno_convergence(
    const std::vector<std::string> &grid_names,
    const std::tuple<double, double> &expected_rate,
    const HybridWENOParams &params) {

  auto f = [](const XYZ &x) {
    //    double g = zisa::sin(2.0 * zisa::pi * (x[0] + x[1] + x[2]));
    auto sigma = 0.2;
    auto r = zisa::norm(x - XYZ{0.5, 0.5, 0.0});
    double g = zisa::exp(-zisa::pow<2>(r / sigma));
    return euler_var_t{g, 0.5 * g, 0.25 * g, 0.4 * g, 0.35 * g};
  };

  std::vector<double> resolution;
  std::vector<double> l1_errors;
  std::vector<double> linf_errors;
  std::vector<std::string> debug_filenames;

  auto quad_deg = zisa::max(2, max_order(params) - 1);

  for (int grid_level = 0; grid_level < grid_names.size(); ++grid_level) {

    auto [grid, sf]
        = load_grid_with_stencils(grid_names, grid_level, params, quad_deg);

    debug_filenames.push_back(string_format(
        "__unit_tests-weno_ao-rate-%d--%d.h5", grid_level, rand()));
    save_grid_with_stencil(debug_filenames.back(), *grid, sf);

    auto n_cells = grid->n_cells;
    auto rc = array<RC, 1>(shape_t<1>{n_cells});
    for (auto i : cell_indices(*grid)) {
      rc[i] = RC(grid, sf[i], i, params);
    }

    auto [l1_err, linf_err] = zisa::compute_errors(*grid, f, rc);

    resolution.push_back(largest_circum_radius(*grid));
    l1_errors.push_back(l1_err);
    linf_errors.push_back(linf_err);
  }

  auto rates = convergence_rates(resolution, l1_errors);

  auto title = string_format("RC = %s", type_name<RC>().c_str());
  auto desc_params = indent_block(1, zisa::to_string(params));
  auto debug_info
      = "debug_files: \n" + indent_block(1, join(debug_filenames, "\n"));

  for (int_t i = 0; i < l1_errors.size(); ++i) {
    auto err_str = string_format(
        "err[%d] = %e, res = %e", i, l1_errors[i], resolution[i]);

    INFO(string_format(
        "%s\n%s\n%s\n%s", title.c_str(), err_str.c_str(), desc_params.c_str()));
    CHECK(zisa::isreal(l1_errors[i]));
  }

  for (int_t i = 0; i < rates.size(); ++i) {
    auto err_str = string_format("err[%d] = %e, rate[%d] = %e, res = %e",
                                 i + 1,
                                 l1_errors[i + 1],
                                 i,
                                 rates[i],
                                 resolution[i + 1]);
    err_str = indent_block(1, err_str);
    INFO(string_format(
        "%s\n%s\n%s", title.c_str(), err_str.c_str(), desc_params.c_str()));
    CHECK(is_inside_interval(rates[i], expected_rate));
  }
}

namespace detail {

inline std::vector<XYZ> make_stability_points_triangular(const Grid &grid,
                                                         int_t i) {
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

inline std::vector<XYZ> make_stability_points_tetrahedral(const Grid &grid,
                                                          int_t i) {
  auto edge_rules = std::vector<TriangularRule>();
  auto volume_rules = std::vector<TetrahedralRule>();

  for (int_t deg = 1; deg <= 4; ++deg) {
    edge_rules.push_back(make_triangular_rule(deg));
    volume_rules.push_back(make_tetrahedral_rule(deg));
  }

  auto points = std::vector<XYZ>();

  auto tet = tetrahedron(grid, i);
  for (int_t k = 0; k < grid.max_neighbours; ++k) {
    auto tri = face(tet, k);

    for (auto &&edge_rule : edge_rules) {
      for (auto &&x : edge_rule.points) {
        points.push_back(coord(tri, x));
      }
    }
  }

  for (auto &&volume_rule : volume_rules) {
    for (auto &&x : volume_rule.points) {
      points.push_back(coord(tet, x));
    }
  }

  return points;
};

inline std::vector<XYZ> make_stability_points(const Grid &grid, int_t i) {
  if (grid.is_triangular()) {
    return make_stability_points_triangular(grid, i);
  } else if (grid.is_tetrahedral()) {
    return make_stability_points_tetrahedral(grid, i);
  }
  LOG_ERR("Implement first.");
}

inline void set_qbar_local(array<euler_var_t, 1> &qbar_local,
                           const HybridWENO &rc,
                           const array<double, 2> &u) {

  const auto &l2g = rc.local2global();

  for (int_t ii = 0; ii < l2g.size(); ++ii) {
    int i = l2g[ii];

    for (int_t k = 0; k < u.shape(1); ++k) {
      qbar_local[ii][k] = u(i, k);
    }
  }
}

}

template <class RC>
void test_hybrid_weno_stability(const std::vector<std::string> &grid_names,
                                const HybridWENOParams &params) {

  using scaling_t = zisa::UnityScaling;
  auto scaling = scaling_t{};

  auto quad_deg = zisa::max(1, max_order(params) - 1);
  double tol = 5e-6;

  for (int grid_level = 0; grid_level < grid_names.size(); ++grid_level) {
    auto grid_name = grid_names[grid_level];

    auto [grid, sf]
        = load_grid_with_stencils(grid_names, grid_level, params, quad_deg);

    auto debug_filename = string_format(
        "__unit_tests-weno_ao-stability-%d--%d.h5", grid_level, rand());
    save_grid_with_stencil(debug_filename, *grid, sf);

    constexpr int_t n_vars = 5;

    auto f = [grid = grid](const XYZ &x) {
      double z = (grid->n_dims() == 2 ? 0.0 : 0.5);
      auto d = zisa::norm(x - XYZ{0.5, 0.5, 0.5});
      return (d < 0.3 ? 10.0 : 1.0);
    };

    auto u = AllVariables({grid->n_cells, int_t(n_vars), int_t(0)});
    for (auto i : cell_indices(*grid)) {
      for (int_t k = 0; k < n_vars; ++k) {
        u.cvars(i, k) = f(grid->cell_centers(i));
      }
    }

    auto n_polys = params.linear_weights.size();
    auto max_stencil_size = 1024;
    using cvars_t = euler_var_t;
    auto qbar_local = array<cvars_t, 1>(shape_t<1>{max_stencil_size});
    auto polys = array<WENOPoly, 1>(shape_t<1>{n_polys});
    auto rhs
        = array<double, 2>(shape_t<2>{max_stencil_size, WENOPoly::n_vars()});

    for (auto i : cell_indices(*grid)) {
      auto stencil_family = sf[i];
      auto rc = RC(grid, stencil_family, i, params);

      detail::set_qbar_local(qbar_local, rc, u.cvars);

      auto weno_poly = rc.reconstruct(rhs, polys, qbar_local);

      auto points = detail::make_stability_points(*grid, i);
      for (auto &&x : points) {

        auto approx = weno_poly(x);
        auto exact = Cartesian<n_vars>(u.cvars(i));

        // clang-format off
        INFO(string_format("[%d] %s - %s = %s\n%s\n%s\n%s\n%s",
                           i,
                           zisa::to_string(approx).c_str(),
                           zisa::to_string(exact).c_str(),
                           zisa::to_string(zisa::Cartesian<n_vars>(approx - exact)).c_str(),
                           string_format("x = %s", format_as_list(x).c_str()).c_str(),
                           string_format("x_center = %s", format_as_list(grid->cell_centers[i]).c_str()).c_str(),
                           string_format("debug_file = %s", debug_filename.c_str()).c_str(),
                           type_name<RC>().c_str(),
                           zisa::to_string(params).c_str()
                           ));
        // clang-format on

        if (!almost_equal(approx, exact, tol)) {
          for (int_t k = 0; k < n_polys; ++k) {
            auto approx = polys(k)(x);
            PRINT(approx);

            const auto &s = stencil_family[k];
            auto i_local = array<double, 1>(shape_t<1>(s.size()));
            auto q_local = array<double, 1>(shape_t<1>(s.size()));
            for (int_t ii = 0; ii < s.size(); ++ii) {
              auto i = s.global(ii);
              q_local[ii] = u.cvars(i, 0);
              i_local[ii] = i;
            }
            PRINT(format_as_list(i_local));
            PRINT(format_as_list(q_local));
          }
        }

        REQUIRE(almost_equal(approx, exact, tol));
      }
    }
  }
}

} // namespace zisa
#endif /* end of include guard */
