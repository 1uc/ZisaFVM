// SPDX-License-Identifier: MIT
// Copyright (c) 2021 ETH Zurich, Luc Grosheintz-Laval

#ifndef EULER_EXPERIMENT_H_0ICHW
#define EULER_EXPERIMENT_H_0ICHW

#include <zisa/config.hpp>
#include <zisa/experiments/numerical_experiment.hpp>
#include <zisa/flux/hllc.hpp>
#include <zisa/io/load_snapshot.hpp>
#include <zisa/math/reference_solution.hpp>
#include <zisa/model/cfl_condition.hpp>
#include <zisa/model/euler_factory.hpp>
#include <zisa/model/local_eos_state.hpp>
#include <zisa/parallelization/halo_exchange.hpp>
#include <zisa/reconstruction/global_reconstruction.hpp>

namespace zisa {

template <class EOS, class Gravity>
class EulerExperiment : public TypicalNumericalExperiment {
private:
  using super = TypicalNumericalExperiment;

protected:
  using eos_t = EOS;
  using gravity_t = Gravity;
  using euler_t = Euler;
  using cvars_t = typename euler_t::cvars_t;
  using flux_t = HLLCBatten<eos_t>;
  using scaling_t = EulerScaling<eos_t>;

public:
  explicit EulerExperiment(const InputParameters &params)
      : super(params),
        euler(make_euler(params)),
        gravity(make_gravity<Gravity>(params)) {}

  EulerExperiment(const InputParameters &params,
                  std::shared_ptr<euler_t> euler);

protected:
  virtual void do_post_run(const std::shared_ptr<AllVariables> &u1) override;
  virtual void do_post_process() override;

  virtual std::shared_ptr<RateOfChange> choose_fvm_rate_of_change();
  virtual std::shared_ptr<RateOfChange> choose_rate_of_change() override;
  virtual std::shared_ptr<RateOfChange> choose_flux_bc() override;
  virtual std::shared_ptr<HaloExchange> choose_halo_exchange();
  virtual std::shared_ptr<Visualization> compute_visualization() override;
  virtual std::shared_ptr<DataSource>
  compute_data_source(std::shared_ptr<FNG> fng);
  virtual std::pair<std::string, std::string> compute_restart_datafile();

  virtual std::shared_ptr<SanityCheck> choose_sanity_check() override;
  virtual std::shared_ptr<CFLCondition> choose_cfl_condition() override;
  virtual AllVariablesDimensions choose_all_variable_dims() override;
  virtual std::pair<std::shared_ptr<AllVariables>,
                    std::shared_ptr<AllVariables>>
  load_initial_conditions() override;

  /// This is used for down-sampling the reference solution.
  virtual std::function<std::shared_ptr<Grid>(const std::string &, int_t)>
  choose_grid_factory();

protected:
  template <class Equilibrium, class RC>
  std::shared_ptr<RateOfChange> choose_physical_rate_of_change();

  template <class Equilibrium, class RC>
  std::shared_ptr<RateOfChange>
  choose_flux_loop(const std::shared_ptr<
                   EulerGlobalReconstruction<Equilibrium, RC, scaling_t>>
                       &global_reconstruction);

  template <class Equilibrium, class RC>
  std::shared_ptr<RateOfChange> choose_gravity_source_loop(
      const std::shared_ptr<
          EulerGlobalReconstruction<Equilibrium, RC, scaling_t>>
          &global_reconstruction);

  template <class Equilibrium, class RC>
  std::shared_ptr<RateOfChange> choose_heating_source_loop(
      const std::shared_ptr<
          EulerGlobalReconstruction<Equilibrium, RC, scaling_t>>
          &global_reconstruction);

  template <class Equilibrium>
  std::shared_ptr<RateOfChange> deduce_reconstruction();

  std::shared_ptr<ReferenceSolution>
  deduce_reference_solution(const AllVariables &u1);

  template <class Equilibrium, class Scaling>
  std::shared_ptr<ReferenceSolution> deduce_reference_solution_eq(
      const AllVariables &u1,
      const std::shared_ptr<
          EulerGlobalReconstruction<Equilibrium, CWENO_AO, Scaling>> &grc);

  template <class Equilibrium, class RC>
  auto choose_reconstruction() -> decltype(auto);

  template <class Equilibrium, class RC, class RCParams>
  auto choose_reconstruction(const RCParams &rc_params) -> decltype(auto);

  HybridWENOParams choose_weno_reference_params() const;
  LocalRCParams choose_local_rc_params() const;

  std::shared_ptr<LocalEOSState<EOS>> choose_local_eos() {
    if (local_eos_ == nullptr) {
      local_eos_ = compute_local_eos();
    }
    return local_eos_;
  }

  virtual std::shared_ptr<LocalEOSState<EOS>> compute_local_eos() {
    auto grid = choose_grid();
    return compute_local_eos(grid->n_cells);
  }

  virtual std::shared_ptr<LocalEOSState<EOS>> compute_local_eos(int_t n_cells) {
    return make_local_eos<EOS>(n_cells, params);
  }

protected:
  std::shared_ptr<euler_t> euler;
  std::shared_ptr<LocalEOSState<eos_t>> local_eos_ = nullptr;
  std::shared_ptr<gravity_t> gravity;
  std::shared_ptr<GlobalReconstruction<euler_var_t>> grc_ = nullptr;
};

} // namespace zisa

#endif /* end of include guard */
