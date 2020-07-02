#ifndef ZISA_STELLAR_CONVECTION_HPP_CGIEYJ
#define ZISA_STELLAR_CONVECTION_HPP_CGIEYJ

#if ZISA_HAS_HELMHOLTZ_EOS == 1
#include <zisa/config.hpp>

#include <zisa/experiments/euler_experiment.hpp>
#include <zisa/math/linear_interpolation.hpp>
#include <zisa/model/helmholtz_eos.hpp>

namespace zisa {

class StellarConvection : public EulerExperiment<HelmholtzEOS, RadialGravity> {
private:
  using super = EulerExperiment<HelmholtzEOS, RadialGravity>;

protected:
  using eos_t = typename super::eos_t;
  using gravity_t = typename super::gravity_t;
  using euler_t = typename super::euler_t;

public:
  StellarConvection(const InputParameters &params) : super(params) {}

protected:
  virtual std::shared_ptr<AllVariables> compute_initial_conditions() override;

  virtual std::function<bool(const Grid &grid, int_t i)>
  boundary_mask() const override;
};

}

#else
#include <zisa/experiments/numerical_experiment.hpp>

namespace zisa {
class StellarConvection : public InvalidNumericalExperiment {
public:
  StellarConvection(const InputParameters &)
      : InvalidNumericalExperiment("Missing `ZISA_HAS_HELMHOLTZ_EOS=1`.") {}
};

}

#endif
#endif // ZISA_STELLAR_CONVECTION_HPP
