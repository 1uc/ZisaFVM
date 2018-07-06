/* Selection of different timestepping schemes.
 *
 * Authors: Luc Grosheintz <forbugrep@zoho.com>
 *    Date: 2015-03-02
 */
#ifndef RUNGE_KUTTA_H_MQEL5YKN
#define RUNGE_KUTTA_H_MQEL5YKN
#include <zisa/config.hpp>

#include <zisa/boundary/boundary_condition.hpp>
#include <zisa/model/all_variables.hpp>
#include <zisa/ode/rate_of_change.hpp>
#include <zisa/ode/time_integration.hpp>

namespace zisa {

/// Rate of change buffers.
class TendencyBuffers {
public:
  TendencyBuffers(const AllVariablesDimensions &dims, int_t n_stages);

  /// Get a reference to the `i`-th buffer.
  AllVariables &operator[](int_t i);

  /// Get a const reference to the `i`-th buffer.
  const AllVariables &operator[](int_t i) const;

  /// Iterator to first buffer.
  std::vector<AllVariables>::iterator begin();

  /// Iterator to past-the-end.
  std::vector<AllVariables>::iterator end();

private:
  std::vector<AllVariables> buffers;
  device_type device = device_type::cpu;
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
 *
 *  Depending on the coefficients the scheme may be explicit or implicit.
 *
 *  NOTE: This tableau only supports explicit schemes, mainly because zisa
 *  only supports explicit time integration and some of the coefficients must
 *  live on the GPU, hence checking them after construction is more
 *  cumbersome.
 */
struct ButcherTableau {
  /// Allocate and initialize the Butcher tableau.
  /** Memory will be allocated on the accelerator as needed, i.e. a and b
   *  but never c. We only consider consistent schemes.
   */
  ButcherTableau(const std::vector<std::vector<double>> &a,
                 const std::vector<double> &b);

public:
  std::vector<array<double, 1>> a; /// `a`-coefficients in accelerator memory.
  array<double, 1> b;              /// `b`-coefficients in accelerator memory.
  array<double, 1> c;              /// `c`-coefficients in host memory

  int_t n_stages; /// number of stages

private:
  void allocate();
  void assign(const std::vector<std::vector<double>> &a,
              const std::vector<double> &b);

private:
  device_type device = device_type::cpu;
};

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
   *  @param cfl_number
   *    The CFL-number controlling stability of the time-integration.
   */
  RungeKutta(const std::shared_ptr<RateOfChange> &rate_of_change,
             const std::shared_ptr<BoundaryCondition> &bc,
             ButcherTableau tableau,
             const AllVariablesDimensions &dims,
             double cfl_number);

  virtual ~RungeKutta() = default;

  virtual std::shared_ptr<AllVariables> compute_step(
      const std::shared_ptr<AllVariables> &u0, double t, double dt) override;

  virtual void boundary_condition(AllVariables &u0, double t) const override;

  virtual double
  pick_time_step(const AllVariables &all_variables) const override;
  virtual double pick_time_step(const AllVariables &all_variables,
                                double dt) const override;

protected:
  std::string assemble_description(int indent, const std::string &detail) const;

protected:
  ButcherTableau tableau;
  std::shared_ptr<RateOfChange> rate_of_change;
  std::shared_ptr<BoundaryCondition> bc;

  std::shared_ptr<AllVariables> ux;
  TendencyBuffers tendency_buffers;

  double cfl_number;

protected:
  device_type device = device_type::cpu;
};

/// Perform the Runge-Kutta step: `u1 = u0 + dt*(a1*k1 + ... + as*ks)`.
/** This interface will dispatch to the requested accelerator.
 */
void runge_kutta_sum(AllVariables &u1,
                     const AllVariables &u0,
                     const TendencyBuffers &tendency_buffers,
                     const array<double, 1> &coeffs,
                     double dt);

/// The simplest first order time integrator.
class ForwardEuler : public RungeKutta {
private:
  using super = RungeKutta;

public:
  ForwardEuler(const std::shared_ptr<RateOfChange> &rate_of_change,
               const std::shared_ptr<BoundaryCondition> &bc,
               const AllVariablesDimensions &dims);

  virtual std::string str(int indent) const override;
};

/// Second order strong stability preserving Runge-Kutta scheme.
class SSP2 : public RungeKutta {
private:
  using super = RungeKutta;

public:
  SSP2(const std::shared_ptr<RateOfChange> &rate_of_change,
       const std::shared_ptr<BoundaryCondition> &bc,
       const AllVariablesDimensions &dims);

  virtual std::string str(int indent) const override;
};

/// Third order strong stability preserving Runge-Kutta scheme.
class SSP3 : public RungeKutta {
private:
  using super = RungeKutta;

public:
  SSP3(const std::shared_ptr<RateOfChange> &rate_of_change,
       const std::shared_ptr<BoundaryCondition> &bc,
       const AllVariablesDimensions &dims);

  virtual std::string str(int indent) const override;
};

/// Third order Wicker-Skamarock time-integrator.
class Wicker : public RungeKutta {
private:
  using super = RungeKutta;

public:
  Wicker(const std::shared_ptr<RateOfChange> &rate_of_change,
         const std::shared_ptr<BoundaryCondition> &bc,
         const AllVariablesDimensions &dims);

  virtual std::string str(int indent) const override;
};

/// The classical Runge-Kutta method.
class RK4 : public RungeKutta {
private:
  using super = RungeKutta;

public:
  RK4(const std::shared_ptr<RateOfChange> &rate_of_change,
      const std::shared_ptr<BoundaryCondition> &bc,
      const AllVariablesDimensions &dims);

  virtual std::string str(int indent) const override;
};

/// The 6-stage fifth-order part of the Fehlberg method.
class Fehlberg : public RungeKutta {
private:
  using super = RungeKutta;

public:
  Fehlberg(const std::shared_ptr<RateOfChange> &rate_of_change,
           const std::shared_ptr<BoundaryCondition> &bc,
           const AllVariablesDimensions &dims);

  virtual std::string str(int indent) const override;
};

} // namespace zisa

// Convince Doxygen to grab my functions as well.
//! @file
#endif /* end of include guard: RUNGE_KUTTA_H_MQEL5YKN */
