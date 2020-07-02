#ifndef EULER_FACTORY_H_1VICG
#define EULER_FACTORY_H_1VICG

#include <type_traits>
#include <utility>

using std::declval;

#include <zisa/cli/input_parameters.hpp>
#include <zisa/model/euler.hpp>
#include <zisa/utils/type_name.hpp>

namespace zisa {

JankaEOSParams make_janka_eos_params(const InputParameters &params);

template <class EOS>
EOS make_eos(const InputParameters & /* params */);

template <>
IdealGasEOS make_eos<IdealGasEOS>(const InputParameters &params);

template <>
JankaEOS make_eos<JankaEOS>(const InputParameters &params);

#if ZISA_HAS_HELMHOLTZ_EOS == 1
template <>
HelmholtzEOS make_eos<HelmholtzEOS>(const InputParameters &params);
#endif

template <class Gravity>
Gravity make_gravity(const InputParameters &);

template <>
NoGravity make_gravity<NoGravity>(const InputParameters &input_params);

template <>
ConstantGravityRadial
make_gravity<ConstantGravityRadial>(const InputParameters &input_params);

template <>
PolytropeGravityRadial
make_gravity<PolytropeGravityRadial>(const InputParameters &input_params);

template <>
PolytropeGravityWithJumpRadial make_gravity<PolytropeGravityWithJumpRadial>(
    const InputParameters &input_params);

template <>
RadialGravity make_gravity<RadialGravity>(const InputParameters &input_params);

template <class Model>
std::shared_ptr<Model> make_euler(const InputParameters &params) {
  using EOS = typename Model::eos_t;
  using Gravity = typename Model::gravity_t;

  auto eos = make_eos<EOS>(params);
  auto gravity = make_gravity<Gravity>(params);

  return std::make_shared<Euler<EOS, Gravity>>(eos, gravity);
}

Euler<IdealGasEOS, ConstantGravityRadial> make_default_euler();

} // namespace zisa
#endif
