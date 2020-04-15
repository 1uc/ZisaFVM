#ifndef NO_EQUILIBRIUM_H_0EJAC
#define NO_EQUILIBRIUM_H_0EJAC

#include <zisa/math/cartesian.hpp>
#include <zisa/model/euler_variables.hpp>

namespace zisa {

struct NoEquilibrium {
  template <class... Args>
  explicit NoEquilibrium(Args &&... /* args */) {}

  RhoE extrapolate(const EnthalpyEntropy &, const XYZ &, const XYZ &) const {
    return {0.0, 0.0};
  }

  std::string str(int /* verbose */ = 0) const { return "NoEquilibrium"; }
};

} // namespace zisa
#endif /* end of include guard */
