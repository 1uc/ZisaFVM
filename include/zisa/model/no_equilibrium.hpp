#ifndef NO_EQUILIBRIUM_H_0EJAC
#define NO_EQUILIBRIUM_H_0EJAC

#include <zisa/model/ideal_gas_eos.hpp>
#include <zisa/model/euler_variables.hpp>
#include <zisa/math/cartesian.hpp>

namespace zisa {

struct NoEquilibrium {
  template<class ...Args>
  NoEquilibrium(Args && ... /* args */) {}

};

template <class Gravity>
RhoE extrapolate(const NoEquilibrium & /* eq */,
                 const IdealGasEOS & /* eos */,
                 const EnthalpyEntropy & /* theta */,
                 const XY & /* xy_ref */,
                 const XY & /* xy */) {

  return {0.0, 0.0};
}

} // namespace zisa

#endif /* end of include guard */
