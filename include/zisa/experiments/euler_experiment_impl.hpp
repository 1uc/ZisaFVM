#ifndef EULER_EXPERIMENT_IMPL_H_VDI33
#define EULER_EXPERIMENT_IMPL_H_VDI33

#include "euler_experiment_decl.hpp"

#include <zisa/boundary/equilibrium_flux_bc.hpp>
#include <zisa/boundary/flux_bc.hpp>
#include <zisa/core/flux_loop.hpp>
#include <zisa/core/gravity_source_loop.hpp>
#include <zisa/experiments/down_sample_reference.hpp>
#include <zisa/io/dump_snapshot.hpp>
#include <zisa/io/euler_plots.hpp>
#include <zisa/io/no_visualization.hpp>
#include <zisa/math/reference_solution.hpp>
#include <zisa/model/isentropic_equilibrium.hpp>
#include <zisa/model/local_cfl_condition.hpp>
#include <zisa/model/no_equilibrium.hpp>
#include <zisa/model/sanity_check_for.hpp>
#include <zisa/reconstruction/cweno_ao.hpp>
#include <zisa/reconstruction/weno_ao.hpp>
#include <zisa/utils/parse_duration.hpp>

namespace zisa {

template <class EOS, class Gravity>
void EulerExperiment<EOS, Gravity>::do_post_run(
    const std::shared_ptr<AllVariables> &u1) {

  if (!has_key(params, "reference")) {
    // Post processing the reference solution is not requested.
    return;
  }

  auto reference_solution = deduce_reference_solution(u1);

  std::vector<std::string> coarse_grid_paths
      = params["reference"]["coarse_grids"];

  down_sample_euler_reference(
      *reference_solution, coarse_grid_paths, "reference.h5");

  HDF5SerialReader reader(file_name_generator->steady_state_filename);
  auto u_delta = std::make_shared<AllVariables>(
      AllVariables::load(reader, all_labels<euler_var_t>()));

  for (int_t i = 0; i < u_delta->size(); ++i) {
    (*u_delta)[i] = (*u1)[i] - (*u_delta)[i];
  }

  auto delta
      = deduce_reference_solution_eq<NoEquilibrium>(u_delta, NoEquilibrium{});
  down_sample_euler_reference(*delta, coarse_grid_paths, "delta.h5");
}

template <class EOS, class Gravity>
void EulerExperiment<EOS, Gravity>::do_post_process() {
  grid = choose_grid();
  file_name_generator = choose_file_name_generator();

  auto data_filename = find_last_data_file(*file_name_generator);
  auto reader = HDF5SerialReader(data_filename);
  auto u1 = std::make_shared<AllVariables>(
      AllVariables::load(reader, all_labels<euler_var_t>()));

  do_post_run(u1);
}

template <class EOS, class Gravity>
template <class Equilibrium>
std::shared_ptr<ReferenceSolution>
EulerExperiment<EOS, Gravity>::deduce_reference_solution_eq(
    const std::shared_ptr<AllVariables> &u1, const Equilibrium &eq) const {

  return std::make_shared<EulerReferenceSolution<Equilibrium>>(grid, u1, eq);
}

template <class EOS, class Gravity>
std::shared_ptr<ReferenceSolution>
EulerExperiment<EOS, Gravity>::deduce_reference_solution(
    const std::shared_ptr<AllVariables> &u1) const {
  if (params["reference"]["equilibrium"] == "constant") {
    return deduce_reference_solution_eq(u1, NoEquilibrium{});
  }

  if (params["reference"]["equilibrium"] == "isentropic") {
    int_t quad_deg = params["quadrature"]["volume"];

    auto eq = IsentropicEquilibrium(euler, quad_deg);
    return deduce_reference_solution_eq(u1, eq);
  }

  LOG_ERR("Failed to deduce reference solutions.");
}

template <class EOS, class Gravity>
std::shared_ptr<CFLCondition>
EulerExperiment<EOS, Gravity>::choose_cfl_condition() {
  LOG_ERR_IF(!has_key(params, "ode"), "Missing section 'ode'.");
  LOG_ERR_IF(!has_key(params["ode"], "cfl_number"),
             "Missing section 'ode/cfl_number'.");

  double cfl_number = params["ode"]["cfl_number"];
  return std::make_shared<LocalCFL<euler_t>>(grid, euler, cfl_number);
}

template <class EOS, class Gravity>
std::shared_ptr<Visualization>
EulerExperiment<EOS, Gravity>::choose_visualization() {
  if (params["io"]["mode"] == "none") {
    return std::make_shared<NoVisualization>();
  } else if (params["io"]["mode"] == "opengl") {
    auto delay = parse_duration_ms(params["io"]["opengl"]["delay"]);
    return std::make_shared<EulerPlots>(*grid, delay);
  } else if (params["io"]["mode"] == "hdf5") {
    const auto &fng = file_name_generator;
    return std::make_shared<SerialDumpSnapshot<euler_t>>(euler, fng);
  }

  LOG_ERR("Implement missing case.");
}

template <class EOS, class Gravity>
std::shared_ptr<SanityCheck>
EulerExperiment<EOS, Gravity>::choose_sanity_check() {
  return std::make_shared<SanityCheckFor<euler_t>>(euler);
}

template <class EOS, class Gravity>
AllVariablesDimensions
EulerExperiment<EOS, Gravity>::choose_all_variable_dims() {
  return {grid->n_cells, euler_t::cvars_t::size(), int_t(0)};
}

template <class EOS, class Gravity>
std::shared_ptr<RateOfChange> EulerExperiment<EOS, Gravity>::choose_flux_bc() {

  std::string flux_bc = params["flux-bc"]["mode"];

  if (flux_bc == "constant") {
    return std::make_shared<FluxBC<euler_t>>(euler, grid);
  }

  if (flux_bc == "isentropic") {
    auto qr = choose_edge_rule();
    using eq_t = IsentropicEquilibrium<eos_t, gravity_t>;
    auto eq = eq_t(euler, params["quadrature"]["volume"]);

    return std::make_shared<EquilibriumFluxBC<eq_t, euler_t>>(
        euler, eq, grid, qr);
  }

  LOG_ERR(
      string_format("Unknown flux boundary condition. [%s]", flux_bc.c_str()));
}

template <class EOS, class Gravity>
template <class Equilibrium, class RC>
std::shared_ptr<RateOfChange>
EulerExperiment<EOS, Gravity>::choose_physical_rate_of_change() {
  auto rc = choose_reconstruction<Equilibrium, RC>();
  auto fvm_change = choose_flux_loop<Equilibrium, RC>(rc);
  auto source_change = choose_gravity_source_loop<Equilibrium, RC>(rc);

  return std::make_shared<SumRatesOfChange>(fvm_change, source_change);
}

template <class EOS, class Gravity>
template <class Equilibrium, class RC>
std::shared_ptr<RateOfChange> EulerExperiment<EOS, Gravity>::choose_flux_loop(
    const std::shared_ptr<EulerGlobalReconstruction<Equilibrium, RC>> &rc) {

  auto edge_rule = choose_edge_rule();
  return std::make_shared<FluxLoop<Equilibrium, RC, euler_t, flux_t>>(
      grid, euler, rc, edge_rule);
}

template <class EOS, class Gravity>
template <class Equilibrium, class RC>
std::shared_ptr<RateOfChange>
EulerExperiment<EOS, Gravity>::choose_gravity_source_loop(
    const std::shared_ptr<EulerGlobalReconstruction<Equilibrium, RC>> &rc) {

  int_t edge_deg = params["quadrature"]["edge"];
  int_t volume_deg = params["quadrature"]["volume"];

  auto edge_rule = choose_edge_rule();
  return std::make_shared<GravitySourceLoop<Equilibrium, RC, euler_t>>(
      grid, euler, rc, edge_deg, volume_deg);
}

template <class EOS, class Gravity>
template <class Equilibrium>
std::shared_ptr<RateOfChange>
EulerExperiment<EOS, Gravity>::deduce_reconstruction() {
  std::string reconstruction = params["reconstruction"]["mode"];

  if (reconstruction == "WENO-AO") {
    return choose_physical_rate_of_change<Equilibrium, WENO_AO>();
  } else if (reconstruction == "CWENO-AO") {
    return choose_physical_rate_of_change<Equilibrium, WENO_AO>();
  }

  LOG_ERR("Failed to deduce reconstruction.");
}

template <class EOS, class Gravity>
std::shared_ptr<RateOfChange>
EulerExperiment<EOS, Gravity>::choose_rate_of_change() {
  auto fvm_change = choose_fvm_rate_of_change();
  return aggregate_rates_of_change({fvm_change});
}

template <class EOS, class Gravity>
std::shared_ptr<RateOfChange>
EulerExperiment<EOS, Gravity>::choose_fvm_rate_of_change() {
  std::string equilibrium = params["well-balancing"]["mode"];

  if (equilibrium == "isentropic") {
    return deduce_reconstruction<IsentropicEquilibrium<eos_t, gravity_t>>();
  } else if (equilibrium == "constant") {
    return deduce_reconstruction<NoEquilibrium>();
  }

  LOG_ERR("Failed to deduce well-balancing.");
}

template <class EOS, class Gravity>
template <class Equilibrium, class RC>
std::shared_ptr<EulerGlobalReconstruction<Equilibrium, RC>>
EulerExperiment<EOS, Gravity>::choose_reconstruction() {
  LOG_ERR_IF(!has_key(params, "reconstruction"),
             "Missing section 'reconstruction'.");
  auto rc_params = params["reconstruction"];

  return choose_reconstruction<Equilibrium, RC>(rc_params);
}

template <class EOS, class Gravity>
template <class Equilibrium, class RC, class RCParams>
std::shared_ptr<EulerGlobalReconstruction<Equilibrium, RC>>
EulerExperiment<EOS, Gravity>::choose_reconstruction(
    const RCParams &rc_params) {
  auto hybrid_weno_params
      = HybridWENOParams(StencilFamilyParams(rc_params["orders"],
                                             rc_params["biases"],
                                             rc_params["overfit_factors"]),
                         rc_params["linear_weights"],
                         rc_params["smoothness_indicator"]["epsilon"],
                         rc_params["smoothness_indicator"]["exponent"]);

  int_t quad_deg = params["quadrature"]["volume"];

  auto eq = Equilibrium{euler, quad_deg};
  return std::make_shared<EulerGlobalReconstruction<Equilibrium, RC>>(
      grid, hybrid_weno_params, eq);
}

} // namespace zisa
#endif /* end of include guard */
