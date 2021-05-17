#ifndef RUNGE_KUTTA_H_MQEL5YKN
#define RUNGE_KUTTA_H_MQEL5YKN

#include <zisa/config.hpp>

#include <functional>

#include <zisa/boundary/boundary_condition.hpp>
#include <zisa/model/all_variables.hpp>
#include <zisa/ode/rate_of_change.hpp>
#include <zisa/ode/time_integration.hpp>

namespace zisa {

/// Rate of change buffers.
class TendencyBuffers {
private:
  std::vector<AllVariables> buffers;

public:
  TendencyBuffers(const AllVariablesDimensions &dims, int_t n_stages);

  /// Get a reference to the `i`-th buffer.
  AllVariables &operator[](int_t i);

  /// Get a const reference to the `i`-th buffer.
  const AllVariables &operator[](int_t i) const;

  /// Iterator to first buffer.
  auto begin() -> decltype(buffers.begin());

  /// Iterator to past-the-end.
  auto end() -> decltype(buffers.end());
};

/// The Butcher tableau for Runge-Kutta time-integrators.
/** A Butcher tableau is a compact representation of any Runge-Kutta method.
 *  The ODE is written as
 *      dy/dt{y}(x, t) = L(t, y).
 *
 *  Runge-Kutta method with s stages can be written as:
 *      y^(j) = y^0 + dt sum_i a_{j,i} k_i     (j = 1...s)
 *      k^(j) = L(t + c_j*dt, y^(j))           (j = 1...s)
 *      c_j = sum_i a_{j,i}                    (j = 1...s)
 *
 *      y^1 = y^0 + dt sum_i b_i k_i.
 *
 *  The coefficients should be such that the scheme is explicit.
 */
struct ButcherTableau {
  /// Allocate and initialize the Butcher tableau.
  /** Memory will be allocated on the accelerator as needed, i.e. a and b
   *  but never c. We only consider consistent schemes.
   */
  ButcherTableau(const std::vector<std::vector<double>> &a,
                 const std::vector<double> &b);

public:
  std::vector<array<double, 1>> a;
  array<double, 1> b;
  array<double, 1> c;

  int_t n_stages;

private:
  void allocate();
  void assign(const std::vector<std::vector<double>> &a,
              const std::vector<double> &b);
};

/// Generate the requested Butcher Tableau.
ButcherTableau make_tableau(const std::string &method);

/// Implementation of Runge-Kutta based on `ButcherTableau`.
class RungeKutta : public TimeIntegration {
private:
  using super = TimeIntegration;

public:
  /** This constructor allocates one tendency buffer for each stage.
   *
   *  @param rate_of_change
   *    Computes the rate of change.
   *  @param boundary_condition
   *    Applies the ghost-cell boundary conditions of the PDE.
   *  @param tableau
   *    Butcher tableau determining the Runge-Kutta scheme.
   *  @param dims
   *    Shape of the rate of change buffers.
   */
  RungeKutta(std::shared_ptr<RateOfChange> rate_of_change,
             std::shared_ptr<BoundaryCondition> bc,
             ButcherTableau tableau,
             const AllVariablesDimensions &dims);

  virtual ~RungeKutta() = default;

  virtual std::shared_ptr<AllVariables> compute_step(
      const std::shared_ptr<AllVariables> &u0, double t, double dt) override;

protected:
  void boundary_condition(AllVariables &u0, double t) const;
  std::string assemble_description(const std::string &detail) const;

protected:
  ButcherTableau tableau;
  std::shared_ptr<RateOfChange> rate_of_change;
  std::shared_ptr<BoundaryCondition> bc;

  std::shared_ptr<AllVariables> ux;
  TendencyBuffers tendency_buffers;
};

template <class X>
class StaticRungeKutta {
private:
  using function_t = std::function<X(double, const X &)>;

public:
  StaticRungeKutta(const function_t f, const std::string &method)
      : f(f), tableau(make_tableau(method)), tendencies(tableau.n_stages) {}

  X operator()(const X &x0, double t, double dt) {
    // stage 0
    tendencies[0] = f(t, x0);

    // stages 1, ..., s
    for (int_t stage = 1; stage < tableau.n_stages; ++stage) {
      double tx = t + tableau.c[stage] * dt;

      auto x = runge_kutta_sum(x0, tableau.a[stage], dt);
      tendencies[stage] = f(tx, x);
    }

    return runge_kutta_sum(x0, tableau.b, dt);
  }

protected:
  X runge_kutta_sum(X x, const array<double, 1> &coeffs, double dt) const {
    int_t n_stages = coeffs.shape(0);
    for (int_t stage = 0; stage < n_stages; ++stage) {
      if (coeffs[stage] != 0.0) {
        x += dt * coeffs[stage] * tendencies[stage];
      }
    }

    return x;
  }

private:
  function_t f;
  ButcherTableau tableau;
  std::vector<X> tendencies;
};

/// The simplest first order time integrator.
class ForwardEuler : public RungeKutta {
private:
  using super = RungeKutta;

public:
  ForwardEuler(const std::shared_ptr<RateOfChange> &rate_of_change,
               const std::shared_ptr<BoundaryCondition> &bc,
               const AllVariablesDimensions &dims);

  virtual std::string str() const override;
};

/// Second order strong stability preserving Runge-Kutta scheme.
class SSP2 : public RungeKutta {
private:
  using super = RungeKutta;

public:
  SSP2(const std::shared_ptr<RateOfChange> &rate_of_change,
       const std::shared_ptr<BoundaryCondition> &bc,
       const AllVariablesDimensions &dims);

  virtual std::string str() const override;
};

/// Third order strong stability preserving Runge-Kutta scheme.
class SSP3 : public RungeKutta {
private:
  using super = RungeKutta;

public:
  SSP3(const std::shared_ptr<RateOfChange> &rate_of_change,
       const std::shared_ptr<BoundaryCondition> &bc,
       const AllVariablesDimensions &dims);

  virtual std::string str() const override;
};

/// Third order Wicker-Skamarock time-integrator.
class Wicker : public RungeKutta {
private:
  using super = RungeKutta;

public:
  Wicker(const std::shared_ptr<RateOfChange> &rate_of_change,
         const std::shared_ptr<BoundaryCondition> &bc,
         const AllVariablesDimensions &dims);

  virtual std::string str() const override;
};

/// The classical Runge-Kutta method.
class RK4 : public RungeKutta {
private:
  using super = RungeKutta;

public:
  RK4(const std::shared_ptr<RateOfChange> &rate_of_change,
      const std::shared_ptr<BoundaryCondition> &bc,
      const AllVariablesDimensions &dims);

  virtual std::string str() const override;
};

/// The 6-stage fifth-order part of the Fehlberg method.
class Fehlberg : public RungeKutta {
private:
  using super = RungeKutta;

public:
  Fehlberg(const std::shared_ptr<RateOfChange> &rate_of_change,
           const std::shared_ptr<BoundaryCondition> &bc,
           const AllVariablesDimensions &dims);

  virtual std::string str() const override;
};

} // namespace zisa

// Convince Doxygen to grab my functions as well.
//! @file
#endif /* end of include guard: RUNGE_KUTTA_H_MQEL5YKN */
