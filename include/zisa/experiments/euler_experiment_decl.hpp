#ifndef EULER_EXPERIMENT_H_0ICHW
#define EULER_EXPERIMENT_H_0ICHW

#include <zisa/config.hpp>
#include <zisa/experiments/numerical_experiment.hpp>
#include <zisa/flux/hllc.hpp>
#include <zisa/model/cfl_condition.hpp>
#include <zisa/model/euler_factory.hpp>
#include <zisa/reconstruction/global_reconstruction.hpp>

namespace zisa {

template <class EOS, class Gravity>
class EulerExperiment : public NumericalExperiment {
private:
  using super = NumericalExperiment;

protected:
  using eos_t = EOS;
  using gravity_t = Gravity;
  using euler_t = Euler<eos_t, gravity_t>;
  using cvars_t = typename euler_t::cvars_t;
  using flux_t = HLLCBatten<euler_t>;

public:
  EulerExperiment(const InputParameters &params)
      : super(params), euler(make_euler<euler_t>(params)) {}

protected:
  virtual std::shared_ptr<RateOfChange> choose_fvm_rate_of_change();
  virtual std::shared_ptr<RateOfChange> choose_rate_of_change() override;
  virtual std::shared_ptr<RateOfChange> choose_flux_bc() override;
  virtual std::shared_ptr<Visualization> choose_visualization() override;
  virtual std::shared_ptr<SanityCheck> choose_sanity_check() override;
  virtual std::shared_ptr<CFLCondition> choose_cfl_condition() override;
  virtual AllVariablesDimensions choose_all_variable_dims() override;

private:
  template <class Equilibrium, class RC>
  std::shared_ptr<RateOfChange> choose_physical_rate_of_change();

  template <class Equilibrium, class RC>
  std::shared_ptr<RateOfChange>
  choose_flux_loop(const std::shared_ptr<GlobalReconstruction<Equilibrium, RC>>
                       &global_reconstruction);

  template <class Equilibrium, class RC>
  std::shared_ptr<RateOfChange> choose_gravity_source_loop(
      const std::shared_ptr<GlobalReconstruction<Equilibrium, RC>>
          &global_reconstruction);

  template <class Equilibrium>
  std::shared_ptr<RateOfChange> deduce_reconstruction();

  template <class Equilibrium, class RC>
  std::shared_ptr<GlobalReconstruction<Equilibrium, RC>>
  choose_reconstruction();

protected:
  const euler_t euler;
};

} // namespace zisa

#endif /* end of include guard */
