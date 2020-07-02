#ifndef ISENTROPIC_EQUILIBRIUM_H_5173V
#define ISENTROPIC_EQUILIBRIUM_H_5173V

#include <zisa/math/triangle.hpp>
#include <zisa/model/equilibrium_values.hpp>
#include <zisa/model/euler.hpp>
#include <zisa/model/euler_variables.hpp>

namespace zisa {

template <class EOS, class Gravity>
struct IsentropicEquilibrium {
private:
  using euler_t = Euler<EOS, Gravity>;

public:
  using equilibrium_values_t = isentropic_equilibrium_values_t<EOS>;

public:
  IsentropicEquilibrium() = default;
  explicit IsentropicEquilibrium(std::shared_ptr<euler_t> euler)
      : euler(std::move(euler)) {}

  RhoE extrapolate(equilibrium_values_t theta_ref,
                   const XYZ &x_ref,
                   const XYZ &x) const {

    const auto &gravity = euler->gravity;
    const auto &eos = euler->eos;

    double phi_ref = gravity.phi(x_ref);
    double phi = gravity.phi(x);

    auto h_ref = theta_ref.h();
    auto s_ref = theta_ref.s();

    double h = h_ref + phi_ref - phi;
    double s = s_ref;

    theta_ref.h() = h;
    theta_ref.s() = s;

    return conserved_variables(eos, theta_ref);
  }

  std::shared_ptr<Euler<EOS, Gravity>> euler;
};

} // namespace zisa
#endif /* end of include guard */
