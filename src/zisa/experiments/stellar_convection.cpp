// SPDX-License-Identifier: MIT
// Copyright (c) 2021 ETH Zurich, Luc Grosheintz-Laval

#include <zisa/experiments/stellar_convection.hpp>

#include <random>

#if ZISA_HAS_HELMHOLTZ_EOS == 1

namespace zisa {
std::pair<std::shared_ptr<AllVariables>, std::shared_ptr<AllVariables>>
StellarConvection::compute_initial_conditions() {
  auto grid = choose_grid();
  auto profile
      = std::string(params["experiment"]["initial_conditions"]["profile"]);

  auto reader = HDF5SerialReader(profile);

  auto points = array<double, 1>::load(reader, "x1");
  auto make_interpolation = [&points, &reader](const std::string &key) {
    return NonUniformLinearInterpolation<double>(
        points, array<double, 1>::load(reader, key));
  };

  auto dims = choose_all_variable_dims();
  auto n_cvars = dims.n_cvars;
  auto n_avars = dims.n_avars;

  reader.open_group("conserved");
  auto cvar_interpolation
      = std::vector<NonUniformLinearInterpolation<double>>();
  cvar_interpolation.reserve(n_cvars);

  auto cvar_keys = std::vector<std::string>{
      "density", "momentum_x1", "momentum_x2", "momentum_x3", "tot._energy"};

  for (auto key : cvar_keys) {
    cvar_interpolation.push_back(make_interpolation(key));
  }
  const auto &rho_interpolation = cvar_interpolation[0];

  reader.switch_group("advected");
  auto avar_keys = std::vector<std::string>();
  for (const auto &element : params["euler"]["eos"]["element_keys"]) {
    avar_keys.push_back(std::string(element));
  }

  auto avar_interpolation
      = std::vector<NonUniformLinearInterpolation<double>>();
  avar_interpolation.reserve(n_avars);

  LOG_ERR_IF(n_avars != avar_keys.size(),
             "Incorrect number of adv. variables.");

  for (auto key : avar_keys) {
    avar_interpolation.push_back(make_interpolation(key));
  }

  reader.switch_group("supplementary");
  auto cs_interpolation = make_interpolation("sound_speed");

  auto all_var_dims = choose_all_variable_dims();
  auto u0 = std::make_shared<AllVariables>(all_var_dims);

  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<double> uniform(-1.0, 1.0);

  for (const auto &[i, cell] : cells(*grid)) {
    auto wrap_interpolation = [](const auto &interpolation) {
      return [&interpolation](const XYZ &x) {
        return interpolation(zisa::norm(x));
      };
    };

    for (int_t k = 0; k < n_cvars; k++) {
      u0->cvars(i, k)
          = average(cell, wrap_interpolation(cvar_interpolation[k]));
    }

    double r = zisa::norm(grid->cell_centers[i]);
    for (int_t k = 1; k < 4; k++) {
      double km = 1e5;
      if (1.2e4 * km < r && r < 2.2e4 * km) {
        double vk = uniform(gen) * 1e-3 * cs_interpolation(r);
        u0->cvars(i, k) += cvar_interpolation[0](r) * vk;
      }
    }

    auto mass_weighted = [&rho_interpolation](const auto &interpolation) {
      return [&rho_interpolation, &interpolation](const XYZ &x) {
        double r = zisa::norm(x);
        return rho_interpolation(r) * interpolation(r);
      };
    };

    for (int_t k = 0; k < n_avars; k++) {
      u0->avars(i, k) = average(cell, mass_weighted(avar_interpolation[k]));
    }
  }

  auto local_eos = choose_local_eos();
  local_eos->compute(*u0);

  auto vis = choose_visualization();
  vis->steady_state(*u0);
  vis->wait();

  return {u0, u0};
}

std::pair<std::shared_ptr<AllVariables>, std::shared_ptr<AllVariables>>
StellarConvection::load_initial_conditions() {
  auto [u0, steady_state] = super::load_initial_conditions();

  auto local_eos = choose_local_eos();
  local_eos->compute(*u0);

  return {u0, steady_state};
}

std::pair<std::shared_ptr<AllVariables>, std::shared_ptr<AllVariables>>
IdealStellarConvection::compute_initial_conditions() {
  auto grid = choose_grid();
  auto profile
      = std::string(params["experiment"]["initial_conditions"]["profile"]);

  auto reader = HDF5SerialReader(profile);

  auto points = array<double, 1>::load(reader, "x1");
  auto make_interpolation = [&points, &reader](const std::string &key) {
    return NonUniformLinearInterpolation<double>(
        points, array<double, 1>::load(reader, key));
  };

  auto dims = choose_all_variable_dims();
  //  auto n_cvars = dims.n_cvars;
  auto n_avars = dims.n_avars;

  reader.open_group("conserved");
  auto rho_interp = make_interpolation("density");

  reader.switch_group("primitive");
  auto p_interp = make_interpolation("pressure");

  reader.switch_group("advected");
  auto avar_keys = std::vector<std::string>();
  for (const auto &element : params["euler"]["eos"]["element_keys"]) {
    avar_keys.push_back(std::string(element));
  }

  auto avar_interpolation
      = std::vector<NonUniformLinearInterpolation<double>>();
  avar_interpolation.reserve(n_avars);

  LOG_ERR_IF(n_avars != avar_keys.size(),
             "Incorrect number of adv. variables.");

  for (auto key : avar_keys) {
    avar_interpolation.push_back(make_interpolation(key));
  }

  reader.switch_group("supplementary");
  auto cs_interpolation = make_interpolation("sound_speed");

  auto all_var_dims = choose_all_variable_dims();
  auto u0 = std::make_shared<AllVariables>(all_var_dims);

  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<double> uniform(-1.0, 1.0);

  for (const auto &[i, cell] : cells(*grid)) {
    auto wrap_interpolation = [](const auto &interpolation) {
      return [&interpolation](const XYZ &x) {
        return interpolation(zisa::norm(x));
      };
    };

    double rho_bar = average(cell, wrap_interpolation(rho_interp));
    double p_bar = average(cell, wrap_interpolation(p_interp));
    double gamma = 5.0 / 3.0;

    u0->cvars(i, 0) = rho_bar;
    u0->cvars(i, 1) = 0.0;
    u0->cvars(i, 2) = 0.0;
    u0->cvars(i, 3) = 0.0;
    u0->cvars(i, 4) = p_bar / (gamma - 1.0);

    double r = zisa::norm(grid->cell_centers[i]);
    for (int_t k = 1; k < 4; k++) {
      double km = 1e5;
      if (1.2e4 * km < r && r < 2.2e4 * km) {
        double vk = uniform(gen) * 1e-3 * cs_interpolation(r);
        u0->cvars(i, k) += rho_bar * vk;
      }
    }
    u0->cvars(i, 4)
        += 0.5
           * (zisa::pow<2>(u0->cvars(i, 1)) + zisa::pow<2>(u0->cvars(i, 2))
              + zisa::pow<2>(u0->cvars(i, 3)))
           / rho_bar;

    auto mass_weighted = [&rho_interp](const auto &interpolation) {
      return [&rho_interp, &interpolation](const XYZ &x) {
        double r = zisa::norm(x);
        return rho_interp(r) * interpolation(r);
      };
    };

    for (int_t k = 0; k < n_avars; k++) {
      u0->avars(i, k) = average(cell, mass_weighted(avar_interpolation[k]));
    }
  }

  auto local_eos = choose_local_eos();
  local_eos->compute(*u0);

  auto vis = choose_visualization();
  vis->steady_state(*u0);
  vis->wait();

  return {u0, u0};
}

std::function<bool(const Grid &grid, int_t i)>
StellarConvection::boundary_mask() const {

  return [](const Grid &grid, int_t i) {
    double km = 1e5 * 1.0;

    double r_inner = 0.0 * km;
    double r_outer = 4e4 * km;
    double r = zisa::norm(grid.cell_centers[i]);
    return (r < r_inner) || (r_outer < r);
  };
}

AllVariablesDimensions StellarConvection::choose_all_variable_dims() {
  auto grid = choose_grid();
  auto n_cells = grid->n_cells;
  auto n_avars = params["euler"]["eos"]["charge_number"].size();

  return {n_cells, 5, n_avars};
}

std::function<bool(const Grid &grid, int_t i)>
IdealStellarConvection::boundary_mask() const {

  return [](const Grid &grid, int_t i) {
    double km = 1e5 * 1.0;

    double r_inner = 5e3 * km;
    double r_outer = 4e4 * km;
    double r = zisa::norm(grid.cell_centers[i]);
    return (r < r_inner) || (r_outer < r);
  };
}

AllVariablesDimensions IdealStellarConvection::choose_all_variable_dims() {
  auto grid = choose_grid();
  auto n_cells = grid->n_cells;
  auto n_avars = params["euler"]["eos"]["charge_number"].size();

  return {n_cells, 5, n_avars};
}

}

#endif
