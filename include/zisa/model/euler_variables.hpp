/*
 *
 * Authors: Luc Grosheintz <forbugrep@zoho.com>
 *    Date: 2017-11-06
 */

#ifndef EULER_VARIABLES_H_4B77E
#define EULER_VARIABLES_H_4B77E

#include <tuple>

#include <zisa/config.hpp>
#include <zisa/math/cartesian.hpp>
#include <zisa/math/edge.hpp>

namespace zisa {

/// Extra variables in Euler equations.
struct euler_xvar_t {
  double p; ///< pressure
  double a; ///< speed of sound
};

/// Conserved variables in euler eq.
/** The conserved variables are (in this order):
 *    - density
 *    - momentum (x, y, z) components.
 *    - total energy
 */
struct euler_var_t : public Cartesian<5> {
  using xvars_t = euler_xvar_t;

  using Cartesian<5>::Cartesian;
  using Cartesian<5>::operator=;

  ANY_DEVICE_INLINE static euler_var_t zeros() {
    return euler_var_t(Cartesian<5>::zeros());
  }

  static std::string labels(int i);
};

ANY_DEVICE_INLINE double kinetic_energy(const euler_var_t &u) {
  return 0.5 * (u[1] * u[1] + u[2] * u[2] + u[3] * u[3]) / u[0];
}

ANY_DEVICE_INLINE double internal_energy(const euler_var_t &u) {
  return u[4] - kinetic_energy(u);
}

void coord_transform(euler_var_t &u, const Edge &edge);
void inv_coord_transform(euler_var_t &u, const Edge &edge);

struct RhoE : public Cartesian<2> {
  using Cartesian<2>::Cartesian;

  ANY_DEVICE_INLINE double rho() const { return (*this)[0]; }
  ANY_DEVICE_INLINE double &rho() { return (*this)[0]; }

  ANY_DEVICE_INLINE double E() const { return (*this)[1]; }
  ANY_DEVICE_INLINE double &E() { return (*this)[1]; }
};

struct RhoT : public Cartesian<2> {
  using Cartesian<2>::Cartesian;

  ANY_DEVICE_INLINE double rho() const { return (*this)[0]; }
  ANY_DEVICE_INLINE double &rho() { return (*this)[0]; }

  ANY_DEVICE_INLINE double T() const { return (*this)[1]; }
  ANY_DEVICE_INLINE double &T() { return (*this)[1]; }
};

struct RhoP : public Cartesian<2> {
  using Cartesian<2>::Cartesian;

  ANY_DEVICE_INLINE double rho() const { return (*this)[0]; }
  ANY_DEVICE_INLINE double &rho() { return (*this)[0]; }

  ANY_DEVICE_INLINE double p() const { return (*this)[1]; }
  ANY_DEVICE_INLINE double &p() { return (*this)[1]; }
};

struct PressureEntropy : public Cartesian<2> {
  using Cartesian<2>::Cartesian;

  ANY_DEVICE_INLINE double p() const { return (*this)[0]; }
  ANY_DEVICE_INLINE double &p() { return (*this)[0]; }

  ANY_DEVICE_INLINE double s() const { return (*this)[1]; }
  ANY_DEVICE_INLINE double &s() { return (*this)[1]; }
};

struct RhoEntropy : public Cartesian<2> {
  using Cartesian<2>::Cartesian;

  ANY_DEVICE_INLINE double rho() const { return (*this)[0]; }
  ANY_DEVICE_INLINE double &rho() { return (*this)[0]; }

  ANY_DEVICE_INLINE double s() const { return (*this)[1]; }
  ANY_DEVICE_INLINE double &s() { return (*this)[1]; }
};

struct EnthalpyEntropy : public Cartesian<2> {
  using Cartesian<2>::Cartesian;

  ANY_DEVICE_INLINE double &h() { return (*this)[0]; }
  ANY_DEVICE_INLINE double h() const { return (*this)[0]; }

  ANY_DEVICE_INLINE double &K() { return (*this)[1]; }
  ANY_DEVICE_INLINE double K() const { return (*this)[1]; }
};

struct Momentum : public Cartesian<3> {
  using Cartesian<3>::Cartesian;
};

inline Momentum momentum(const euler_var_t &u) {
  return Momentum{u[1], u[2], u[3]};
}

} // namespace zisa

#define ZISA_ENABLE_STRUCTRED_BINDINGS(ClassName)                              \
  namespace std {                                                              \
  template <size_t i>                                                          \
  struct tuple_element<i, zisa::ClassName> {                                   \
    using type = double;                                                       \
  };                                                                           \
                                                                               \
  template <>                                                                  \
  struct tuple_size<zisa::ClassName>                                           \
      : public integral_constant<size_t, zisa::ClassName::size()> {};          \
  }                                                                            \
                                                                               \
  namespace zisa {                                                             \
  template <int i>                                                             \
  auto get(const ClassName &cls) {                                             \
    return cls[i];                                                             \
  }                                                                            \
  }

ZISA_ENABLE_STRUCTRED_BINDINGS(euler_var_t);
ZISA_ENABLE_STRUCTRED_BINDINGS(RhoE);
ZISA_ENABLE_STRUCTRED_BINDINGS(RhoT);
ZISA_ENABLE_STRUCTRED_BINDINGS(RhoP);
ZISA_ENABLE_STRUCTRED_BINDINGS(RhoEntropy);
ZISA_ENABLE_STRUCTRED_BINDINGS(PressureEntropy);
ZISA_ENABLE_STRUCTRED_BINDINGS(EnthalpyEntropy);
ZISA_ENABLE_STRUCTRED_BINDINGS(Momentum);

#undef ZISA_ENABLE_STRUCTRED_BINDINGS
#endif /* end of include guard */
