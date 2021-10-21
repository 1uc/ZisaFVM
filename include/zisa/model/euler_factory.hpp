// SPDX-License-Identifier: MIT
// Copyright (c) 2021 ETH Zurich, Luc Grosheintz-Laval

#ifndef EULER_FACTORY_H_1VICG
#define EULER_FACTORY_H_1VICG

#include <type_traits>
#include <utility>

#include <zisa/cli/input_parameters.hpp>
#include <zisa/model/euler.hpp>
#include <zisa/model/local_eos_state.hpp>
#include <zisa/utils/type_name.hpp>

namespace zisa {

JankaEOSParams make_janka_eos_params(const InputParameters &params);

template <class EOS>
std::shared_ptr<LocalEOSState<EOS>>
make_local_eos(int_t n_cells, const InputParameters & /* params */);

template <>
std::shared_ptr<LocalEOSState<IdealGasEOS>>
make_local_eos<IdealGasEOS>(int_t n_cells, const InputParameters &params);

#if ZISA_HAS_HELMHOLTZ_EOS == 1
template <>
std::shared_ptr<LocalEOSState<HelmholtzEOS>>
make_local_eos<HelmholtzEOS>(int_t n_cells, const InputParameters &params);
#endif

template <class Gravity>
std::shared_ptr<Gravity> make_gravity(const InputParameters &);

template <>
std::shared_ptr<NoGravity>
make_gravity<NoGravity>(const InputParameters &input_params);

template <>
std::shared_ptr<ConstantGravityRadial>
make_gravity<ConstantGravityRadial>(const InputParameters &input_params);

template <>
std::shared_ptr<PolytropeGravityRadial>
make_gravity<PolytropeGravityRadial>(const InputParameters &input_params);

template <>
std::shared_ptr<PolytropeGravityWithJumpRadial>
make_gravity<PolytropeGravityWithJumpRadial>(
    const InputParameters &input_params);

template <>
std::shared_ptr<RadialGravity>
make_gravity<RadialGravity>(const InputParameters &input_params);

std::shared_ptr<Euler> make_euler(const InputParameters &params);

std::tuple<Euler, IdealGasEOS, ConstantGravityRadial> make_default_euler();

} // namespace zisa
#endif
