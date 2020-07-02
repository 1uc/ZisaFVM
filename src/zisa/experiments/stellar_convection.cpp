#include <zisa/experiments/stellar_convection.hpp>

#if ZISA_HAS_HELMHOLTZ_EOS == 1

namespace zisa {

std::shared_ptr<AllVariables> StellarConvection::compute_initial_conditions() {
  auto grid = choose_grid();
  auto profile
      = std::string(params["experiment"]["initial_conditions"]["profile"]);

  auto reader = HDF5SerialReader(profile);

  auto points = array<double, 1>::load(reader, "x1");
  auto make_interpolation = [&points, &reader](const std::string &key) {
    return NonUniformLinearInterpolation<double>(
        points, array<double, 1>::load(reader, key));
  };

  reader.open_group("conserved");
  auto density = make_interpolation("density");
  auto momentum_x1 = make_interpolation("momentum_x1");
  auto momentum_x2 = make_interpolation("momentum_x2");
  auto momentum_x3 = make_interpolation("momentum_x3");
  auto tot_energy = make_interpolation("tot._energy");

  reader.switch_group("auxiliary");
  auto grav_pot = make_interpolation("grav._potential");

  auto all_var_dims = choose_all_variable_dims();
  auto u0 = std::make_shared<AllVariables>(all_var_dims);

  for (const auto &[i, cell] : cells(*grid)) {
    auto wrap_interpolation = [](const auto &interpolation) {
      return [&interpolation](const XYZ &x) {
        return interpolation(zisa::norm(x));
      };
    };

    u0->cvars(i, 0) = average(cell, wrap_interpolation(density));
    u0->cvars(i, 1) = average(cell, wrap_interpolation(momentum_x1));
    u0->cvars(i, 2) = average(cell, wrap_interpolation(momentum_x2));
    u0->cvars(i, 3) = average(cell, wrap_interpolation(momentum_x3));
    u0->cvars(i, 4) = average(cell, wrap_interpolation(tot_energy));
  }

  auto vis = choose_visualization();
  vis->steady_state(*u0);
  vis->wait();

  return u0;
}

std::function<bool(const Grid &grid, int_t i)>
StellarConvection::boundary_mask() const {

  return [](const Grid &grid, int_t i) {
    double km = 1e5 * 1.0;

    double r_inner = 5e3 * km;
    double r_outer = 4e4 * km;
    double r = zisa::norm(grid.cell_centers[i]);
    return (r < r_inner) || (r_outer < r);
  };
}

}

#endif
