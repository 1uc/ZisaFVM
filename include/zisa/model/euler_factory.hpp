#ifndef EULER_FACTORY_H_1VICG
#define EULER_FACTORY_H_1VICG

#include <type_traits>
#include <utility>

using std::declval;

#include <zisa/cli/input_parameters.hpp>
#include <zisa/model/euler.hpp>
#include <zisa/utils/type_name.hpp>

namespace zisa {

template <class EOS>
EOS make_eos(const InputParameters & /* params */) {
  LOG_ERR(string_format("Need to implement this case. [%s]",
                        type_name<EOS>().c_str()));
}

template <>
inline IdealGasEOS make_eos<IdealGasEOS>(const InputParameters &params) {
  assert(has_key(params, "euler"));
  auto euler_params = params["euler"];

  assert(has_key(euler_params, "eos"));
  auto eos_params = euler_params["eos"];

  assert(has_key(eos_params, "gamma"));
  double gamma = eos_params["gamma"];

  assert(has_key(eos_params, "specific-gas-constant"));
  double R = eos_params["specific-gas-constant"];

  return IdealGasEOS(gamma, R);
}

template <class Gravity>
Gravity make_gravity(const InputParameters &) {
  LOG_ERR(string_format("Need to implement this case. [%s]",
                        type_name<Gravity>().c_str()));
}

template <>
ConstantGravityRadial
make_gravity<ConstantGravityRadial>(const InputParameters &input_params);

template <>
PolytropeGravityRadial
make_gravity<PolytropeGravityRadial>(const InputParameters &input_params);

template <class Model>
Model make_euler(const InputParameters &params) {
  using EOS = typename Model::eos_t;
  using Gravity = typename Model::gravity_t;

  auto eos = make_eos<EOS>(params);
  auto gravity = make_gravity<Gravity>(params);

  return Euler<EOS, Gravity>(eos, gravity);
}

Euler<IdealGasEOS, ConstantGravityRadial> make_default_euler();

} // namespace zisa
#endif
