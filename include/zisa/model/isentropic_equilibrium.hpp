#ifndef ISENTROPIC_EQUILIBRIUM_H_5173V
#define ISENTROPIC_EQUILIBRIUM_H_5173V

#include <zisa/math/triangle.hpp>
#include <zisa/model/euler.hpp>
#include <zisa/model/euler_variables.hpp>

namespace zisa {

template <class EOS, class Gravity>
struct IsentropicEquilibrium {
  IsentropicEquilibrium() = default;
  IsentropicEquilibrium(const Euler<EOS, Gravity> &euler, int_t quad_deg)
      : eos(euler.eos), gravity(euler.gravity), quad_deg(quad_deg) {}

  IsentropicEquilibrium(const EOS &eos, const Gravity &gravity, int_t quad_deg)
      : eos(eos), gravity(gravity), quad_deg(quad_deg) {}

  RhoE extrapolate(const EnthalpyEntropy &theta,
                   const XYZ &x_ref,
                   const XYZ &x) const {

    double phi_ref = gravity.phi(x_ref);
    double phi = gravity.phi(x);

    double h = theta.h() + phi_ref - phi;
    double K = theta.s();

    return eos.rhoE(EnthalpyEntropy{h, K});
  }

  EOS eos;
  Gravity gravity;
  int_t quad_deg; // FIXME this is wrong / unneeded.
};

} // namespace zisa
#endif /* end of include guard */
