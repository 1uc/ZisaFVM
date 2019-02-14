#ifndef NO_EQUILIBRIUM_H_0EJAC
#define NO_EQUILIBRIUM_H_0EJAC

#include <zisa/math/cartesian.hpp>
#include <zisa/model/euler_variables.hpp>
#include <zisa/model/ideal_gas_eos.hpp>

namespace zisa {

struct NoEquilibrium {
  template <class... Args>
  NoEquilibrium(Args &&... /* args */) {}
};

template <class Gravity>
RhoE extrapolate(const NoEquilibrium & /* eq */,
                 const IdealGasEOS & /* eos */,
                 const EnthalpyEntropy & /* theta */,
                 const XYZ & /* xy_ref */,
                 const XYZ & /* xy */) {

  return {0.0, 0.0};
}

} // namespace zisa

#endif /* end of include guard */
