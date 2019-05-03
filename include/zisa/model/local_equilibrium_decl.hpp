#ifndef LOCAL_EQUILIBRIUM_DECL_H_Z7M0R
#define LOCAL_EQUILIBRIUM_DECL_H_Z7M0R

#include <zisa/math/cartesian.hpp>
#include <zisa/math/triangle.hpp>
#include <zisa/model/euler_variables.hpp>
#include <zisa/model/isentropic_equilibrium.hpp>
#include <zisa/model/janka_eos.hpp>
#include <zisa/model/no_equilibrium.hpp>

namespace zisa {

template <class Equilibrium>
class LocalEquilibriumBase {
public:
  LocalEquilibriumBase() = default;
  explicit LocalEquilibriumBase(const Equilibrium &equilibrium);
  LocalEquilibriumBase(const Equilibrium &equilibrium,
                       const EnthalpyEntropy &theta_ref,
                       const XYZ &x);

  void solve(const RhoE &rhoE_bar, const Triangle &tri_ref);
  RhoE extrapolate(const XYZ &xy) const;
  RhoE extrapolate(const Triangle &tri) const;

protected:
  EnthalpyEntropy theta = EnthalpyEntropy{};
  XYZ x_ref = XYZ{};
  bool found_equilibrium = false;

  Equilibrium equilibrium;
};

template <class Equilibrium>
class LocalEquilibrium : public LocalEquilibriumBase<Equilibrium> {
private:
  using super = LocalEquilibriumBase<Equilibrium>;

public:
  LocalEquilibrium() = default;
  explicit LocalEquilibrium(const Equilibrium &equilibrium)
      : super(equilibrium) {}
  LocalEquilibrium(const Equilibrium &equilibrium,
                   const EnthalpyEntropy &theta_ref,
                   const XYZ &x)
      : super(equilibrium, theta_ref, x) {}
};

template <>
class LocalEquilibrium<IsentropicEquilibrium<JankaEOS, RadialGravity>>
    : public LocalEquilibriumBase<
          IsentropicEquilibrium<JankaEOS, RadialGravity>> {
private:
  using eq_t = IsentropicEquilibrium<JankaEOS, RadialGravity>;
  using super = LocalEquilibriumBase<eq_t>;

public:
  LocalEquilibrium() = default;
  explicit LocalEquilibrium(const eq_t &equilibrium) : super(equilibrium) {}
  LocalEquilibrium(const eq_t &equilibrium,
                   const EnthalpyEntropy &theta_ref,
                   const XYZ &x)
      : super(equilibrium, theta_ref, x) {}

  void solve(const RhoE &rhoE_bar, const Triangle &tri_ref) {
    const auto &eos = this->equilibrium.euler->eos;
    auto [rho, E] = rhoE_bar;

    E = zisa::max(eos.polytropic_energy(rho), E);

    super::solve(RhoE{rho, E}, tri_ref);
  }
};

template <>
class LocalEquilibrium<NoEquilibrium> {
public:
  LocalEquilibrium() = default;

  template <class... Args>
  explicit LocalEquilibrium(Args &&... /* args */) {}

  template <class... Args>
  inline void solve(Args &&... /* args */) {
    // do nothing
  }

  template <class... Args>
  RhoE extrapolate(Args &&... /* args */) const {
    return {0.0, 0.0};
  }
};

} // namespace zisa
#endif /* end of include guard */
