#ifndef POLYTROPE_H_5K155
#define POLYTROPE_H_5K155

#include <zisa/config.hpp>
#include <zisa/experiments/euler_experiment.hpp>
#include <zisa/model/gravity.hpp>
#include <zisa/model/ideal_gas_eos.hpp>

namespace zisa {

class Polytrope : public EulerExperiment<IdealGasEOS, PolytropeGravityRadial> {
private:
  using super = EulerExperiment<IdealGasEOS, PolytropeGravityRadial>;

protected:
  using eos_t = typename super::eos_t;
  using gravity_t = typename super::gravity_t;
  using euler_t = typename super::euler_t;
  using cvars_t = typename super::cvars_t;

public:
  using super::super;

protected:
  virtual std::shared_ptr<AllVariables> choose_initial_conditions() override;
};

} // namespace zisa

#endif
