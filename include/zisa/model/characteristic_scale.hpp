#ifndef ZISA_CHARACTERISTIC_SCALE_HPP
#define ZISA_CHARACTERISTIC_SCALE_HPP

#include <zisa/config.hpp>

#include <zisa/model/euler_variables.hpp>

namespace zisa {
template <class EULER>
class EulerScaling {
private:
  using euler_t = EULER;
  using cvars_t = euler_var_t;

public:
  EulerScaling() = default;
  explicit EulerScaling(std::shared_ptr<euler_t> euler)
      : euler(std::move(euler)) {}

  cvars_t operator()(const RhoE &rhoE) const {
    auto [rho, E] = rhoE;
    auto xvars = euler->eos.xvars(rhoE);
    double cs = xvars.a;

    LOG_ERR_IF(rho <= 0.0, "Invalid density.");
    LOG_ERR_IF(E <= 0.0, "Invalid energy.");

    return {rho, cs, cs, cs, E};
  }

private:
  std::shared_ptr<euler_t> euler;
};

class UnityScaling {
  using cvars_t = euler_var_t;

public:
  cvars_t operator()(const RhoE &) const { return {1.0, 1.0, 1.0, 1.0, 1.0}; }
};
}

#endif // ZISA_CHARACTERISTIC_SCALE_HPP
