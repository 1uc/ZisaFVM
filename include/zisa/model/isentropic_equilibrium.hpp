#ifndef ISENTROPIC_EQUILIBRIUM_H_5173V
#define ISENTROPIC_EQUILIBRIUM_H_5173V

#include "no_equilibrium.hpp"
#include <zisa/math/triangle.hpp>
#include <zisa/model/equilibrium_values.hpp>
#include <zisa/model/euler.hpp>
#include <zisa/model/euler_variables.hpp>

namespace zisa {

template <class EOS, class Gravity>
struct IsentropicEquilibrium {
public:
  using equilibrium_values_t = isentropic_equilibrium_values_t<EOS>;

public:
  IsentropicEquilibrium() = default;
  explicit IsentropicEquilibrium(std::shared_ptr<EOS> eos,
                                 std::shared_ptr<Gravity> gravity)
      : eos(std::move(eos)), gravity(std::move(gravity)) {}

  RhoE extrapolate(equilibrium_values_t theta_ref,
                   const XYZ &x_ref,
                   const XYZ &x) const {

    double phi_ref = gravity->phi(x_ref);
    double phi = gravity->phi(x);

    auto h_ref = theta_ref.h();
    auto s_ref = theta_ref.s();

    double h = h_ref + phi_ref - phi;
    double s = s_ref;

    theta_ref.h() = h;
    theta_ref.s() = s;

    return conserved_variables(*eos, theta_ref);
  }

  std::shared_ptr<EOS> eos;
  std::shared_ptr<Gravity> gravity;
};

namespace detail {

template <class Equilibrium>
class MakeEquilibrium;

template <class EOS, class Gravity>
class MakeEquilibrium<IsentropicEquilibrium<EOS, Gravity>> {
public:
  static IsentropicEquilibrium<EOS, Gravity>
  make_obj(const std::shared_ptr<EOS> &eos,
           const std::shared_ptr<Gravity> &gravity) {

    return IsentropicEquilibrium(eos, gravity);
  }
};

template <>
class MakeEquilibrium<NoEquilibrium> {
public:
  template <class... T>
  static NoEquilibrium make_obj(T &&...) {
    return NoEquilibrium{};
  }
};

}

template <class Equilibrium, class EOS, class Gravity>
Equilibrium make_equilibrium(const std::shared_ptr<EOS> &eos,
                             const std::shared_ptr<Gravity> &gravity) {

  return detail::MakeEquilibrium<Equilibrium>::make_obj(eos, gravity);
}

} // namespace zisa
#endif /* end of include guard */
