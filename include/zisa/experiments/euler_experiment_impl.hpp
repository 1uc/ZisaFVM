#ifndef EULER_EXPERIMENT_IMPL_H_VDI33
#define EULER_EXPERIMENT_IMPL_H_VDI33

#include "euler_experiment_decl.hpp"
#include <zisa/boundary/flux_bc.hpp>
#include <zisa/core/flux_loop.hpp>
#include <zisa/io/euler_plots.hpp>
#include <zisa/io/no_visualization.hpp>
#include <zisa/model/isentropic_equilibrium.hpp>
#include <zisa/model/local_cfl_condition.hpp>
#include <zisa/model/no_equilibrium.hpp>
#include <zisa/model/sanity_check_for.hpp>
#include <zisa/reconstruction/cweno_ao.hpp>
#include <zisa/reconstruction/weno_ao.hpp>
#include <zisa/utils/parse_duration.hpp>

namespace zisa {

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
  return std::make_shared<FluxBC<euler_t>>(euler, grid);
}

template <class EOS, class Gravity>
template <class Equilibrium>
std::shared_ptr<RateOfChange>
EulerExperiment<EOS, Gravity>::deduce_reconstruction() {
  auto edge_rule = choose_edge_rule();

  std::string reconstruction = params["reconstruction"]["mode"];

  if (reconstruction == "WENO-AO") {
    auto rc = choose_reconstruction<Equilibrium, WENO_AO>();
    return std::make_shared<FluxLoop<Equilibrium, WENO_AO, euler_t, flux_t>>(
        grid, euler, rc, edge_rule);
  } else if (reconstruction == "CWENO-AO") {
    auto rc = choose_reconstruction<Equilibrium, CWENO_AO>();
    return std::make_shared<FluxLoop<Equilibrium, CWENO_AO, euler_t, flux_t>>(
        grid, euler, rc, edge_rule);
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
  } else if (equilibrium == "naive") {
    return deduce_reconstruction<NoEquilibrium>();
  }

  LOG_ERR("Failed to deduce well-balancing.");
}

template <class EOS, class Gravity>
template <class Equilibrium, class RC>
std::shared_ptr<GlobalReconstruction<Equilibrium, RC>>
EulerExperiment<EOS, Gravity>::choose_reconstruction() {
  LOG_ERR_IF(!has_key(params, "reconstruction"),
             "Missing section 'reconstruction'.");
  auto rc_params = params["reconstruction"];

  auto hybrid_weno_params
      = HybridWENO_Params(StencilFamilyParams(rc_params["orders"],
                                              rc_params["biases"],
                                              rc_params["overfit_factors"]),
                          rc_params["linear_weights"],
                          rc_params["smoothness_indicator"]["epsilon"],
                          rc_params["smoothness_indicator"]["exponent"]);

  int_t quad_deg = params["quadrature"]["volume"];

  auto eq = Equilibrium{euler.eos, euler.gravity, quad_deg};
  return std::make_shared<GlobalReconstruction<Equilibrium, RC>>(
      grid, hybrid_weno_params, eq);
}

} // namespace zisa
#endif /* end of include guard */
