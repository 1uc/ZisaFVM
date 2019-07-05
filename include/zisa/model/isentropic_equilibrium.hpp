#ifndef ISENTROPIC_EQUILIBRIUM_H_5173V
#define ISENTROPIC_EQUILIBRIUM_H_5173V

#include <zisa/math/triangle.hpp>
#include <zisa/model/euler.hpp>
#include <zisa/model/euler_variables.hpp>

namespace zisa {

template <class EOS, class Gravity>
struct IsentropicEquilibrium {
private:
  using euler_t = Euler<EOS, Gravity>;

public:
  IsentropicEquilibrium() = default;
  explicit IsentropicEquilibrium(std::shared_ptr<euler_t> euler)
      : euler(std::move(euler)) {}

  RhoE extrapolate(const EnthalpyEntropy &theta,
                   const XYZ &x_ref,
                   const XYZ &x) const {

    const auto &gravity = euler->gravity;
    const auto &eos = euler->eos;

    double phi_ref = gravity.phi(x_ref);
    double phi = gravity.phi(x);

    double h = theta.h() + phi_ref - phi;
    double K = theta.s();

    return eos.rhoE(EnthalpyEntropy{h, K});
  }

  std::shared_ptr<Euler<EOS, Gravity>> euler;
  //  int_t quad_deg; // FIXME this is wrong / unneeded.
};

} // namespace zisa
#endif /* end of include guard */
