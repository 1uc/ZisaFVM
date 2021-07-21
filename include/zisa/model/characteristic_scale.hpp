// SPDX-License-Identifier: MIT
// Copyright (c) 2021 ETH Zurich, Luc Grosheintz-Laval

#ifndef ZISA_CHARACTERISTIC_SCALE_HPP
#define ZISA_CHARACTERISTIC_SCALE_HPP

#include <zisa/config.hpp>

#include <zisa/model/euler_variables.hpp>
#include <zisa/model/local_eos_state.hpp>

namespace zisa {
template <class EOS>
class EulerScaling {
private:
  using eos_t = EOS;
  using cvars_t = euler_var_t;

public:
  EulerScaling() = default;
  explicit EulerScaling(std::shared_ptr<EOS> eos) : eos(std::move(eos)) {}

  cvars_t operator()(const RhoE &rhoE) const {
    auto [rho, E] = rhoE;
    auto xvars = eos->xvars(rhoE);
    double cs = xvars.a;

    LOG_ERR_IF(rho <= 0.0, string_format("Invalid density. [%e]", rho));
    LOG_ERR_IF(E <= 0.0, string_format("Invalid energy. [%e]", E));

    return {rho, cs, cs, cs, E};
  }

private:
  std::shared_ptr<eos_t> eos;
};

class UnityScaling {
  using cvars_t = euler_var_t;

public:
  template <class... T>
  explicit UnityScaling(T &&...) {}

  cvars_t operator()(const RhoE &) const { return {1.0, 1.0, 1.0, 1.0, 1.0}; }
};

template <class EOS>
class LocalEulerScaling {
public:
  LocalEulerScaling(std::shared_ptr<LocalEOSState<EOS>> local_eos)
      : local_eos(std::move(local_eos)) {}

  EulerScaling<EOS> operator()(int_t i) const {
    return EulerScaling<EOS>((*local_eos)(i));
  }

private:
  std::shared_ptr<LocalEOSState<EOS>> local_eos;
};

class LocalUnityScaling {
public:
  UnityScaling operator()(int_t /* i */) const { return UnityScaling{}; }
};

}

#endif // ZISA_CHARACTERISTIC_SCALE_HPP
