// SPDX-License-Identifier: MIT
// Copyright (c) 2021 ETH Zurich, Luc Grosheintz-Laval

#ifndef ZISA_EQUILIBRIUM_VALUES_HPP_BIUWQ
#define ZISA_EQUILIBRIUM_VALUES_HPP_BIUWQ

#include <zisa/model/ideal_gas_eos.hpp>

namespace zisa {

template <class EOS>
struct IsentropicEquilibriumValuesTraits;

template <>
struct IsentropicEquilibriumValuesTraits<IdealGasEOS> {
  using equilibrium_values_t = EnthalpyEntropy;
};

template <class EOS>
using isentropic_equilibrium_values_t =
    typename IsentropicEquilibriumValuesTraits<EOS>::equilibrium_values_t;

ANY_DEVICE_INLINE isentropic_equilibrium_values_t<IdealGasEOS>
isentropic_equilibrium_values(const IdealGasEOS &,
                              const EnthalpyEntropy &theta,
                              const RhoT &) {
  return theta;
}

ANY_DEVICE_INLINE RhoE conserved_variables(const IdealGasEOS &eos,
                                           const EnthalpyEntropy &theta) {
  return eos.rhoE(theta);
}

ANY_DEVICE_INLINE std::pair<RhoE, euler_xvar_t>
full_variables(const IdealGasEOS &eos, const EnthalpyEntropy &theta) {
  auto rhoE = eos.rhoE(theta);
  return {rhoE, eos.xvars(rhoE)};
}

}

#if ZISA_HAS_HELMHOLTZ_EOS == 1
#include <zisa/model/euler_variables.hpp>
#include <zisa/model/helmholtz_eos.hpp>

namespace zisa {

struct HelmholtzIsentropicEquilibriumValues {
  EnthalpyEntropy enthalpy_entropy;
  RhoT rhoT_guess;

  double &h() { return enthalpy_entropy.h(); }
  double h() const { return enthalpy_entropy.h(); }

  double &s() { return enthalpy_entropy.s(); }
  double s() const { return enthalpy_entropy.s(); }
};

template <>
struct IsentropicEquilibriumValuesTraits<HelmholtzEOS> {
  using equilibrium_values_t = HelmholtzIsentropicEquilibriumValues;
};

ANY_DEVICE_INLINE
isentropic_equilibrium_values_t<HelmholtzEOS>
isentropic_equilibrium_values(const HelmholtzEOS &,
                              const EnthalpyEntropy &theta,
                              const RhoT &rhoT_guess) {
  return {theta, rhoT_guess};
}

ANY_DEVICE_INLINE
RhoE conserved_variables(const HelmholtzEOS &eos,
                         const HelmholtzIsentropicEquilibriumValues &theta) {
  return eos.rhoE(theta.enthalpy_entropy, theta.rhoT_guess);
}

ANY_DEVICE_INLINE std::pair<RhoE, euler_xvar_t>
full_variables(const HelmholtzEOS &eos,
               const HelmholtzIsentropicEquilibriumValues &theta) {

  auto full_xvars
      = eos.full_extra_variables(theta.enthalpy_entropy, theta.rhoT_guess);
  return {eos.rhoE(full_xvars), eos.xvars(full_xvars)};
}

}

#endif
#endif // ZISA_EQUILIBRIUM_VALUES_HPP
