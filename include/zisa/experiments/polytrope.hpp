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
  virtual std::shared_ptr<AllVariables> compute_initial_conditions() override;
  std::shared_ptr<AllVariables> compute_initial_conditions(double amp,
                                                           double width);

  std::shared_ptr<BoundaryCondition> compute_boundary_condition() override;

  /// This is used for down-sampling the reference solution.
  std::function<std::shared_ptr<Grid>(const std::string &, int_t)>
  choose_grid_factory() override;

  void enforce_cell_flags(Grid &grid) const override;
};

class JankaBump : public EulerExperiment<IdealGasEOS, PolytropeGravityRadial> {
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
  virtual std::shared_ptr<AllVariables> compute_initial_conditions() override;
  std::shared_ptr<AllVariables> compute_initial_conditions(double amp,
                                                           double width);
};

} // namespace zisa

#endif
