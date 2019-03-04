#ifndef EULER_FACTORY_H_BDRPI
#define EULER_FACTORY_H_BDRPI

#include <zisa/model/euler_factory.hpp>

namespace zisa {

template <>
ConstantGravityRadial
make_gravity<ConstantGravityRadial>(const InputParameters &params) {
  LOG_ERR_IF(params["euler"]["gravity"]["mode"] != "constant",
             "Incompatible gravity.");

  return ConstantGravityRadial(double(params["euler"]["gravity"]["g"]));
}

template <>
PolytropeGravityRadial
make_gravity<PolytropeGravityRadial>(const InputParameters &params) {
  LOG_ERR_IF(params["euler"]["gravity"]["mode"] != "polytrope",
             "Incompatible gravity.");

  return PolytropeGravityRadial();
}

Euler<IdealGasEOS, ConstantGravityRadial> make_default_euler() {
  return Euler{IdealGasEOS{1.2, 2.3}, ConstantGravityRadial{0.99}};
}

} // namespace zisa

#endif /* end of include guard */
