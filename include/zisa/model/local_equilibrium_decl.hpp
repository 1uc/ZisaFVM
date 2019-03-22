#ifndef LOCAL_EQUILIBRIUM_DECL_H_Z7M0R
#define LOCAL_EQUILIBRIUM_DECL_H_Z7M0R

#include <zisa/math/cartesian.hpp>
#include <zisa/math/triangle.hpp>
#include <zisa/model/euler_variables.hpp>
#include <zisa/model/no_equilibrium.hpp>

namespace zisa {

template <class Equilibrium>
class LocalEquilibrium {
public:
  LocalEquilibrium() = default;
  LocalEquilibrium(const Equilibrium &equilibrium);

  void solve(const RhoE &rhoE_bar, const Triangle &tri_ref);
  RhoE extrapolate(const XYZ &xy) const;
  RhoE extrapolate(const Triangle &tri) const;

private:
  EnthalpyEntropy theta = EnthalpyEntropy{};
  XYZ x_ref = XYZ{};
  bool found_equilibrium = false;

  Equilibrium equilibrium;
};

template <>
class LocalEquilibrium<NoEquilibrium> {
public:
  LocalEquilibrium() = default;

  template <class... Args>
  LocalEquilibrium(Args &&... /* args */) {}

  template <class... Args>
  inline void solve(Args &&... /* args */) {
    return;
  }

  template <class... Args>
  RhoE extrapolate(Args &&... /* args */) const {
    return {0.0, 0.0};
  }
};

} // namespace zisa
#endif /* end of include guard */
