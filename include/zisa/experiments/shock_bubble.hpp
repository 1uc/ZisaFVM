#ifndef SHOCK_BUBBLE_H_RZ1IZ
#define SHOCK_BUBBLE_H_RZ1IZ

#include <zisa/config.hpp>
#include <zisa/experiments/numerical_experiment.hpp>
#include <zisa/flux/hllc.hpp>
#include <zisa/model/cfl_condition.hpp>
#include <zisa/model/euler_factory.hpp>
#include <zisa/model/isentropic_equilibrium.hpp>
#include <zisa/reconstruction/cweno_ao.hpp>
#include <zisa/reconstruction/global_reconstruction.hpp>

namespace zisa {

class ShockBubble : public NumericalExperiment {
private:
  using super = NumericalExperiment;
  using eos_t = IdealGasEOS;
  using gravity_t = ConstantGravityRadial;
  using euler_t = Euler<eos_t, gravity_t>;
  using cvars_t = typename euler_t::cvars_t;
  using flux_t = HLLCBatten<euler_t>;
  // using eq_t = IsentropicEquilibrium<eos_t, gravity_t> ;
  using eq_t = NoEquilibrium;
  using lrc_t = CWENO_AO;
  using grc_t = GlobalReconstruction<eq_t, lrc_t>;

public:
  ShockBubble(const InputParameters &params)
      : super(params), euler(make_euler<euler_t>(params)) {}

protected:
  virtual std::shared_ptr<grc_t> choose_reconstruction();

  virtual std::shared_ptr<AllVariables> choose_initial_conditions() override;
  virtual std::shared_ptr<RateOfChange> choose_rate_of_change() override;
  virtual std::shared_ptr<Visualization> choose_visualization() override;
  virtual std::shared_ptr<SanityCheck> choose_sanity_check() override;
  virtual std::shared_ptr<CFLCondition> choose_cfl_condition() override;
  virtual AllVariablesDimensions choose_all_variable_dims() override;

private:
  euler_t euler;
};

} // namespace zisa

#endif
