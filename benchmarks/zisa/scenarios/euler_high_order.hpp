#ifndef EULER_HIGH_ORDER_H_41DV8
#define EULER_HIGH_ORDER_H_41DV8

#include <zisa/boundary/flux_bc.hpp>
#include <zisa/boundary/no_boundary_condition.hpp>
#include <zisa/core/flux_loop.hpp>
#include <zisa/core/gravity_source_loop.hpp>
#include <zisa/flux/hllc.hpp>
#include <zisa/math/edge_rule.hpp>
#include <zisa/model/euler.hpp>
#include <zisa/ode/runge_kutta.hpp>
#include <zisa/reconstruction/cweno_ao.hpp>
#include <zisa/reconstruction/global_reconstruction.hpp>

namespace zisa {
namespace scenarios {
namespace euler_high_order {

namespace types {

using euler_t = Euler<IdealGasEOS, PolytropeGravityRadial>;
using eq_t = NoEquilibrium;
using rc_t = CWENO_AO;
using global_reconstruction_t = GlobalReconstruction<eq_t, rc_t>;
using flux_loop_t = FluxLoop<eq_t, rc_t, euler_t, HLLCBatten<euler_t>>;
using source_loop_t = GravitySourceLoop<eq_t, rc_t, euler_t>;

} // namespace types

inline std::shared_ptr<Grid> load_grid() {
  return load_gmsh("grids/convergence/unit_square_2.msh");
}

inline types::euler_t make_model() {
  return {IdealGasEOS(/* gamma = */ 1.2, /* R = */ 1.1), PolytropeGravityRadial()};
}

inline types::eq_t make_equilibrium() {
  auto euler = make_model();
  return types::eq_t(euler.eos, euler.gravity);
}

inline std::shared_ptr<types::global_reconstruction_t>
make_global_reconstruction(std::shared_ptr<Grid> &grid) {

  int_t n_vars = 5;
  double eps = 1e-12;
  double s = 4.0;

  auto params = HybridWENO_Params{
      {{4, 2, 2, 2}, {"c", "b", "b", "b"}, {2.0, 1.5, 1.5, 1.5}},
      {10.0, 1.0, 1.0, 1.0},
      eps,
      s};

  return std::make_shared<types::global_reconstruction_t>(
      grid, params, make_equilibrium());
}

inline std::shared_ptr<types::flux_loop_t> make_flux_loop(
    std::shared_ptr<Grid> grid,
    std::shared_ptr<types::global_reconstruction_t> global_reconstruction) {

  auto model = make_model();
  auto edge_rule = EdgeRule(4);

  return std::make_shared<types::flux_loop_t>(
      grid, model, global_reconstruction, edge_rule);
}

inline std::shared_ptr<types::flux_loop_t>
make_flux_loop(std::shared_ptr<Grid> grid) {
  auto global_reconstruction = make_global_reconstruction(grid);
  return make_flux_loop(grid, global_reconstruction);
}

inline std::shared_ptr<types::source_loop_t> make_source_loop(
    std::shared_ptr<Grid> grid,
    std::shared_ptr<types::global_reconstruction_t> global_reconstruction) {

  auto model = make_model();
  int_t edge_deg = 4;
  int_t volume_deg = 4;

  return std::make_shared<types::source_loop_t>(
      grid, model, global_reconstruction, edge_deg, volume_deg);
}

inline std::shared_ptr<types::source_loop_t>
make_source_loop(std::shared_ptr<Grid> grid) {
  auto global_reconstruction = make_global_reconstruction(grid);
  return make_source_loop(grid, global_reconstruction);
}

inline std::shared_ptr<AllVariables> make_all_variables(int_t n_cells) {
  return std::make_shared<AllVariables>(
      AllVariablesDimensions{n_cells, int_t(5), int_t(0)});
}

inline std::shared_ptr<FluxBC<types::euler_t>>
make_flux_bc(const std::shared_ptr<Grid> &grid) {
  return std::make_shared<FluxBC<types::euler_t>>(make_model(), grid);
}

inline std::shared_ptr<ZeroRateOfChange> make_zero_rate_of_change() {
  return std::make_shared<ZeroRateOfChange>();
}

inline std::shared_ptr<TimeIntegration>
make_time_integration(const std::shared_ptr<Grid> &grid) {
  auto zero_change = make_zero_rate_of_change();
  auto fvm_change = make_flux_loop(grid);
  auto flux_bc = make_flux_bc(grid);

  auto roc = std::make_shared<SumRatesOfChange>();
  roc->add_term(zero_change);
  roc->add_term(fvm_change);
  roc->add_term(flux_bc);

  auto no_bc = std::make_shared<NoBoundaryCondition>();
  auto dims = AllVariablesDimensions{grid->n_cells, int_t(5), int_t(0)};

  return std::make_shared<SSP3>(roc, no_bc, dims);
}

inline std::shared_ptr<AllVariables>
make_valid_initial_conditions(const Grid &grid) {
  auto ic = make_all_variables(grid.n_cells);

  for (auto &&[i, tri] : triangles(grid)) {
    auto xy = barycenter(tri);
    auto x = xy[0];
    auto y = xy[1];

    auto rho = 1.0 + 0.1 * x * y;
    auto vx = sin(x) * sin(y);
    auto vy = pow<2>(sin(3.0 * x * y)) * cos(10.0 * y * y);
    auto E = 0.5 * rho * (vx * vx + vy * vy) + 1.0;

    ic->cvars(i) = euler_var_t{rho, rho * vx, rho * vy, 0.0, E};
  }

  return ic;
}

} // namespace euler_high_order
} // namespace scenarios
} // namespace zisa

#endif /* end of include guard */
