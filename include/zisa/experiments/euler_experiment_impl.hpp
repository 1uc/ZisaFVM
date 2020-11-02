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
#include <zisa/model/heating.hpp>
#include <zisa/model/isentropic_equilibrium.hpp>
#include <zisa/model/local_cfl_condition.hpp>
#include <zisa/model/no_equilibrium.hpp>
#include <zisa/model/sanity_check_for.hpp>
#include <zisa/reconstruction/cweno_ao.hpp>
#include <zisa/reconstruction/weno_ao.hpp>
#include <zisa/utils/parse_duration.hpp>

namespace zisa {

template <class EOS, class Gravity>
EulerExperiment<EOS, Gravity>::EulerExperiment(const InputParameters &params,
                                               std::shared_ptr<euler_t> euler_)
    : super(params), euler(std::move(euler_)) {

  if (is_restart()) {
    auto reader = HDF5SerialReader(params["restart"]["file"]);
    LOG_ERR("load EOS & Gravity.");
  }
}

template <class EOS, class Gravity>
void EulerExperiment<EOS, Gravity>::do_post_run(
    const std::shared_ptr<AllVariables> &u1) {

  if (!has_key(params, "reference")) {
    // Post processing the reference solution is not requested.
    return;
  }

  auto reference_solution = deduce_reference_solution(*u1);

  std::vector<std::string> coarse_grid_paths
      = params["reference"]["coarse_grids"];

  auto grid_factory = choose_grid_factory();

  down_sample_euler_reference(
      *reference_solution, coarse_grid_paths, grid_factory, "reference.h5");

  auto fng = choose_file_name_generator();
  auto steady_state_filename = fng->steady_state_filename;
  auto u_delta = std::make_shared<AllVariables>(load_serial<AllVariables>(
      steady_state_filename, all_labels<euler_var_t>()));

  for (int_t i = 0; i < u_delta->size(); ++i) {
    (*u_delta)[i] = (*u1)[i] - (*u_delta)[i];
  }

  auto grid = choose_grid();
  auto local_eos = choose_local_eos();
  auto weno_params = choose_weno_reference_params();
  auto local_rc_params = choose_local_rc_params();
  auto rc = make_reconstruction_array<NoEquilibrium,
                                      CWENO_AO,
                                      UnityScaling,
                                      EOS,
                                      Gravity>(
      grid, weno_params, *local_eos, gravity, local_rc_params);
  auto grc = std::make_shared<
      EulerGlobalReconstruction<NoEquilibrium, CWENO_AO, UnityScaling>>(
      weno_params, std::move(rc));

  auto delta = deduce_reference_solution_eq(*u_delta, grc);
  down_sample_euler_reference(
      *delta, coarse_grid_paths, grid_factory, "delta.h5");
}

template <class EOS, class Gravity>
LocalRCParams EulerExperiment<EOS, Gravity>::choose_local_rc_params() const {
  auto steps_per_recompute
      = int_t(params["reconstruction"]["steps_per_recompute"]);

  auto recompute_threshold
      = double(params["reconstruction"]["recompute_threshold"]);

  return LocalRCParams{steps_per_recompute, recompute_threshold};
}

template <class EOS, class Gravity>
void EulerExperiment<EOS, Gravity>::do_post_process() {
  auto fng = choose_file_name_generator();

  auto data_filename = find_last_data_file(*fng);
  auto reader = HDF5SerialReader(data_filename);
  auto u1 = std::make_shared<AllVariables>(
      AllVariables::load(reader, all_labels<euler_var_t>()));

  do_post_run(u1);
}

template <class EOS, class Gravity>
std::shared_ptr<AllVariables>
EulerExperiment<EOS, Gravity>::load_initial_conditions() {
  std::string datafile = params["restart"]["file"];

  auto reader = HDF5SerialReader(datafile);
  return std::make_shared<AllVariables>(
      AllVariables::load(reader, all_labels<euler_var_t>()));
}

template <class EOS, class Gravity>
template <class Equilibrium, class Scaling>
std::shared_ptr<ReferenceSolution>
EulerExperiment<EOS, Gravity>::deduce_reference_solution_eq(
    const AllVariables &u1,
    const std::shared_ptr<
        EulerGlobalReconstruction<Equilibrium, CWENO_AO, Scaling>> &grc) {

  auto grid = choose_grid();
  return std::make_shared<EulerReferenceSolution<Equilibrium, Scaling>>(
      grid, u1, grc);
}

template <class EOS, class Gravity>
HybridWENOParams
EulerExperiment<EOS, Gravity>::choose_weno_reference_params() const {

  auto grid = choose_grid();
  auto n_dims = grid->n_dims();

  if (n_dims == 2) {
    return HybridWENOParams{
        {{5, 2, 2, 2}, {"c", "b", "b", "b"}, {2.0, 1.5, 1.5, 1.5}},
        {100.0, 1.0, 1.0, 1.0},
        1e-10,
        4.0};
  }

  if (n_dims == 3) {
    return HybridWENOParams{
        {{4, 2, 2, 2, 2}, {"c", "b", "b", "b", "b"}, {3.0, 2.0, 2.0, 2.0, 2.0}},
        {100.0, 1.0, 1.0, 1.0, 1.0},
        1e-10,
        4.0};
  }

  LOG_ERR("Implement first.");
}

template <class EOS, class Gravity>
std::shared_ptr<ReferenceSolution>
EulerExperiment<EOS, Gravity>::deduce_reference_solution(
    const AllVariables &u1) {
  auto local_eos = choose_local_eos();
  auto grid = choose_grid();
  auto weno_params = choose_weno_reference_params();
  auto local_rc_params = choose_local_rc_params();

  if (params["reference"]["equilibrium"] == "constant") {
    LOG_WARN("No, this needs to be local");
    auto rc = make_reconstruction_array<NoEquilibrium,
                                        CWENO_AO,
                                        EulerScaling<eos_t>,
                                        eos_t,
                                        gravity_t>(
        grid, weno_params, *local_eos, gravity, local_rc_params);

    auto grc = std::make_shared<EulerGlobalReconstruction<NoEquilibrium,
                                                          CWENO_AO,
                                                          EulerScaling<eos_t>>>(
        weno_params, std::move(rc));
    return deduce_reference_solution_eq(u1, std::move(grc));
  }

  if (params["reference"]["equilibrium"] == "isentropic") {
    auto rc = make_reconstruction_array<IsentropicEquilibrium<eos_t, gravity_t>,
                                        CWENO_AO,
                                        EulerScaling<eos_t>,
                                        eos_t,
                                        gravity_t>(
        grid, weno_params, *local_eos, gravity, local_rc_params);
    auto grc = std::make_shared<
        EulerGlobalReconstruction<IsentropicEquilibrium<eos_t, gravity_t>,
                                  CWENO_AO,
                                  EulerScaling<eos_t>>>(weno_params,
                                                        std::move(rc));

    return deduce_reference_solution_eq(u1, grc);
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
  auto grid = choose_grid();
  auto local_eos = choose_local_eos();
  return std::make_shared<LocalCFL<eos_t>>(grid, euler, local_eos, cfl_number);
}

template <class EOS, class Gravity>
std::shared_ptr<Visualization>
EulerExperiment<EOS, Gravity>::compute_visualization() {
  auto grid = choose_grid();

  if (params["io"]["mode"] == "none") {
    return std::make_shared<NoVisualization>();
  } else if (params["io"]["mode"] == "opengl") {
    auto delay = parse_duration_ms(params["io"]["opengl"]["delay"]);
    return std::make_shared<EulerPlots>(*grid, delay);
  } else if (params["io"]["mode"] == "hdf5") {
    const auto &fng = choose_file_name_generator();
    auto local_eos = compute_local_eos();
    return std::make_shared<SerialDumpSnapshot<eos_t>>(local_eos, fng);
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
  auto grid = choose_grid();
  return {grid->n_cells, euler_t::cvars_t::size(), choose_n_avars()};
}

template <class EOS, class Gravity>
std::shared_ptr<RateOfChange> EulerExperiment<EOS, Gravity>::choose_flux_bc() {
  std::string flux_bc = params["flux-bc"]["mode"];
  auto grid = choose_grid();
  auto local_eos = choose_local_eos();

  if (flux_bc == "constant") {
    return std::make_shared<FluxBC<eos_t>>(euler, local_eos, grid);
  }

  if (flux_bc == "isentropic") {
    auto qr = choose_edge_rule();
    using eq_t = IsentropicEquilibrium<eos_t, gravity_t>;

    return std::make_shared<EquilibriumFluxBC<eq_t, EOS, Gravity>>(
        euler, local_eos, gravity, grid, qr);
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
  auto heating_change = choose_heating_source_loop<Equilibrium, RC>(rc);

  return std::make_shared<SumRatesOfChange>(
      fvm_change, source_change, heating_change);
}

template <class EOS, class Gravity>
template <class Equilibrium, class RC>
std::shared_ptr<RateOfChange> EulerExperiment<EOS, Gravity>::choose_flux_loop(
    const std::shared_ptr<EulerGlobalReconstruction<Equilibrium, RC, scaling_t>>
        &rc) {
  using grc_t = EulerGlobalReconstruction<Equilibrium, RC, scaling_t>;
  auto grid = choose_grid();
  auto local_eos = choose_local_eos();

  auto edge_rule = choose_edge_rule();
  return std::make_shared<
      FluxLoop<euler_t, flux_t, LocalEOSState<eos_t>, grc_t>>(
      grid, euler, local_eos, rc, edge_rule);
}

template <class EOS, class Gravity>
template <class Equilibrium, class RC>
std::shared_ptr<RateOfChange>
EulerExperiment<EOS, Gravity>::choose_gravity_source_loop(
    const std::shared_ptr<EulerGlobalReconstruction<Equilibrium, RC, scaling_t>>
        &rc) {
  auto grid = choose_grid();
  auto local_eos = choose_local_eos();
  return std::make_shared<GravitySourceLoop<Equilibrium,
                                            RC,
                                            LocalEOSState<eos_t>,
                                            gravity_t,
                                            scaling_t>>(
      grid, local_eos, gravity, rc);
}

template <class EOS, class Gravity>
template <class Equilibrium, class RC>
std::shared_ptr<RateOfChange>
EulerExperiment<EOS, Gravity>::choose_heating_source_loop(
    const std::shared_ptr<EulerGlobalReconstruction<Equilibrium, RC, scaling_t>>
        &rc) {
  auto grid = choose_grid();
  return make_heating_source(grid, rc, params);
}

template <class EOS, class Gravity>
template <class Equilibrium>
std::shared_ptr<RateOfChange>
EulerExperiment<EOS, Gravity>::deduce_reconstruction() {
  std::string reconstruction = params["reconstruction"]["mode"];

  if (reconstruction == "WENO-AO") {
    return choose_physical_rate_of_change<Equilibrium, WENO_AO>();
  } else if (reconstruction == "CWENO-AO") {
    return choose_physical_rate_of_change<Equilibrium, CWENO_AO>();
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
auto EulerExperiment<EOS, Gravity>::choose_reconstruction() -> decltype(auto) {
  LOG_ERR_IF(!has_key(params, "reconstruction"),
             "Missing section 'reconstruction'.");
  auto rc_params = params["reconstruction"];

  return choose_reconstruction<Equilibrium, RC>(rc_params);
}

template <class EOS, class Gravity>
template <class Equilibrium, class RC, class RCParams>
auto EulerExperiment<EOS, Gravity>::choose_reconstruction(
    const RCParams &rc_params) -> decltype(auto) {
  auto hybrid_weno_params
      = HybridWENOParams(StencilFamilyParams(rc_params["orders"],
                                             rc_params["biases"],
                                             rc_params["overfit_factors"]),
                         rc_params["linear_weights"],
                         rc_params["smoothness_indicator"]["epsilon"],
                         rc_params["smoothness_indicator"]["exponent"]);

  auto grid = choose_grid();
  auto stencils = choose_stencils();
  auto local_eos = choose_local_eos();
  auto local_rc_params = choose_local_rc_params();

  auto rc = make_reconstruction_array<Equilibrium,
                                      RC,
                                      EulerScaling<EOS>,
                                      EOS,
                                      Gravity>(grid,
                                               *stencils,
                                               hybrid_weno_params,
                                               *local_eos,
                                               gravity,
                                               local_rc_params);

  return std::make_shared<
      EulerGlobalReconstruction<Equilibrium, RC, scaling_t>>(hybrid_weno_params,
                                                             rc);
}

template <class EOS, class Gravity>
std::function<std::shared_ptr<Grid>(const std::string &, int_t)>
EulerExperiment<EOS, Gravity>::choose_grid_factory() {
  return [](const std::string &grid_name, int_t quad_deg) {
    return load_grid(grid_name, QRDegrees{quad_deg, quad_deg, quad_deg});
  };
}

} // namespace zisa
#endif /* end of include guard */
