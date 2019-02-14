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
  LocalEquilibrium(const Equilibrium &equilibrium, const Triangle &tri_ref);

  void solve(const RhoE &rhoE_bar);
  RhoE extrapolate(const XYZ &xy) const;
  RhoE extrapolate(const Triangle &tri) const;

private:
  EnthalpyEntropy theta;
  bool found_equilibrium;

  Triangle tri_ref;
  Equilibrium equilibrium;
};

template <>
class LocalEquilibrium<NoEquilibrium> {
public:
  LocalEquilibrium() = default;

  template <class... Args>
  LocalEquilibrium(Args &&... /* args */) {}

  inline void solve(const RhoE &/* rhoE_bar */) { return; }
  inline RhoE extrapolate(const XYZ & /* xy */) const { return {0.0, 0.0}; }
  inline RhoE extrapolate(const Triangle &/* tri */) const { return {0.0, 0.0}; }
};

} // namespace zisa
#endif /* end of include guard */
