#ifndef SHOCK_BUBBLE_H_RZ1IZ
#define SHOCK_BUBBLE_H_RZ1IZ

#include <zisa/config.hpp>
#include <zisa/experiments/numerical_experiment.hpp>
#include <zisa/model/euler_factory.hpp>
#include <zisa/flux/hllc.hpp>
#include <zisa/model/cfl_condition.hpp>

namespace zisa {

class ShockBubble : public NumericalExperiment {
private:
  using super = NumericalExperiment;
  using euler_t = Euler<IdealGasEOS, ConstantGravityRadial>;
  using cvars_t = typename euler_t::cvars_t;
  using flux_t = HLLCBatten<euler_t>;

public:
  ShockBubble(const InputParameters &params) : super(params) {
    euler = make_euler<euler_t>(params);
  }

protected:
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
