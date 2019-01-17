#ifndef LOCAL_EQUILIBRIUM_DECL_H_Z7M0R
#define LOCAL_EQUILIBRIUM_DECL_H_Z7M0R

#include <zisa/model/euler_variables.hpp>
#include <zisa/math/cartesian.hpp>
#include <zisa/math/triangle.hpp>

namespace zisa {

template <class Equilibrium>
class LocalEquilibrium {
public:
  LocalEquilibrium(const Equilibrium &equilibrium, const Triangle &tri_ref);

  void solve(const RhoE &rhoE_bar);
  RhoE extrapolate(const XY &xy) const;
  RhoE extrapolate(const Triangle &tri) const;

private:
  EnthalpyEntropy theta;
  bool found_equilibrium;

  Triangle tri_ref;
  Equilibrium equilibrium;
};

} // namespace zisa
#endif /* end of include guard */
