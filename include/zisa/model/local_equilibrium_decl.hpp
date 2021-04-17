#ifndef LOCAL_EQUILIBRIUM_DECL_H_Z7M0R
#define LOCAL_EQUILIBRIUM_DECL_H_Z7M0R

#include <zisa/math/cartesian.hpp>
#include <zisa/math/cell.hpp>
#include <zisa/model/euler_variables.hpp>
#include <zisa/model/isentropic_equilibrium.hpp>
#include <zisa/model/janka_eos.hpp>
#include <zisa/model/no_equilibrium.hpp>

namespace zisa {

template <class Equilibrium>
class LocalEquilibriumBase {
private:
  using xvars_t = euler_xvar_t;

public:
  using equilibrium_values_t = typename Equilibrium::equilibrium_values_t;

public:
  LocalEquilibriumBase() = default;
  explicit LocalEquilibriumBase(const Equilibrium &equilibrium);
  LocalEquilibriumBase(const Equilibrium &equilibrium,
                       const equilibrium_values_t &theta_ref,
                       const XYZ &x);

  void solve(const RhoE &rhoE_bar, const Cell &cell_ref);
  RhoE extrapolate(const XYZ &xy) const;
  std::pair<RhoE, xvars_t> extrapolate_full(const XYZ &xy) const;

  RhoE extrapolate(const Cell &cell) const;

  std::string str(int verbose = 0) const;

protected:
  void solve_exact(const RhoE &rhoE_bar, const Cell &cell_ref);
  void solve_lsq(const RhoE &rhoE_bar, const Cell &cell_ref);

  equilibrium_values_t theta = equilibrium_values_t{};
  XYZ x_ref = XYZ{};
  bool found_equilibrium = false;

  Equilibrium equilibrium;
};

template <class Equilibrium>
class LocalEquilibrium : public LocalEquilibriumBase<Equilibrium> {
private:
  using super = LocalEquilibriumBase<Equilibrium>;

public:
  using equilibrium_values_t = typename super::equilibrium_values_t;

public:
  LocalEquilibrium() = default;
  explicit LocalEquilibrium(const Equilibrium &equilibrium)
      : super(equilibrium) {}
  LocalEquilibrium(const Equilibrium &equilibrium,
                   const equilibrium_values_t &theta_ref,
                   const XYZ &x)
      : super(equilibrium, theta_ref, x) {}
};

template <>
class LocalEquilibrium<NoEquilibrium> {
public:
  LocalEquilibrium() = default;

  template <class... Args>
  explicit LocalEquilibrium(Args &&.../* args */) {}

  template <class... Args>
  inline void solve(Args &&.../* args */) {
    // do nothing
  }

  template <class... Args>
  RhoE extrapolate(Args &&.../* args */) const {
    return {0.0, 0.0};
  }

  template <class... Args>
  std::pair<RhoE, euler_xvar_t> extrapolate_full(Args &&.../* args */) const {
    return {RhoE{0.0, 0.0}, euler_xvar_t{0.0, 0.0}};
  }

  std::string str(int /* verbose */ = 0) const {
    return "LocalEquilibrium<NoEquilibrium>";
  }
};

} // namespace zisa
#endif /* end of include guard */
