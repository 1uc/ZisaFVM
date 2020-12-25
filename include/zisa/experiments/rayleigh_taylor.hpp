#ifndef ZISA_RAYLEIGH_TAYLOR_28379_HPP
#define ZISA_RAYLEIGH_TAYLOR_28379_HPP

#include <zisa/experiments/euler_experiment.hpp>

namespace zisa {

class RayleighTaylor
    : public EulerExperiment<IdealGasEOS, PolytropeGravityRadial> {
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
  virtual std::pair<std::shared_ptr<AllVariables>,
                    std::shared_ptr<AllVariables>>
  compute_initial_conditions() override;
  std::shared_ptr<AllVariables> compute_initial_conditions(double amp,
                                                           double amp_noise,
                                                           double width,
                                                           int n_bumps);

  int_t choose_n_avars() override;

  std::function<bool(const Grid &grid, int_t i)> boundary_mask() const override;
};

}

#endif // ZISA_RAYLEIGH_TAYLOR_HPP
