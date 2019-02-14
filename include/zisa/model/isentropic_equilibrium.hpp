#ifndef ISENTROPIC_EQUILIBRIUM_H_5173V
#define ISENTROPIC_EQUILIBRIUM_H_5173V

#include <zisa/math/triangle.hpp>
#include <zisa/model/euler_variables.hpp>
#include <zisa/model/ideal_gas_eos.hpp>

namespace zisa {

template <class EOS, class Gravity>
struct IsentropicEquilibrium {
  IsentropicEquilibrium() = default;
  IsentropicEquilibrium(const EOS &eos, const Gravity &gravity, int_t quad_deg)
      : eos(eos), gravity(gravity), quad_deg(quad_deg) {}

  EOS eos;
  Gravity gravity;
  int_t quad_deg; // FIXME this is wrong / unneeded.
};

template <class Gravity>
RhoE extrapolate(const IsentropicEquilibrium<IdealGasEOS, Gravity> &eq,
                 const IdealGasEOS &eos,
                 const EnthalpyEntropy &theta,
                 const XYZ &xy_ref,
                 const XYZ &xy) {
  const auto &gravity = eq.gravity;

  double phi_ref = gravity.phi(xy_ref);
  double phi = gravity.phi(xy);

  double h = theta.h() + phi_ref - phi;
  double K = theta.K();

  return eos.rhoE(EnthalpyEntropy{h, K});
}

} // namespace zisa
#endif /* end of include guard */
