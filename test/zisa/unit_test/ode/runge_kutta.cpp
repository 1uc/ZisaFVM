#include <numeric>

#include <zisa/testing/testing_framework.hpp>

#include <zisa/math/basic_functions.hpp>
#include <zisa/model/all_variables.hpp>
#include <zisa/ode/time_integration_factory.hpp>

namespace zisa {

class DoNothingBC : public BoundaryCondition {
public:
  virtual ~DoNothingBC() = default;

  virtual void apply(AllVariables &, double) override { return; }

  virtual std::string str() const override { return "Do nothing mock BC."; }
};

class MockRHS {
public:
  virtual double dudt(double t, double u) const = 0;
  virtual double solution(double t, double u) const = 0;
};

class ConstantRHS : public MockRHS {
public:
  virtual double dudt(double /* t */, double /* u */) const override {
    return 0.1;
  }
  virtual double solution(double t, double u0) const override {
    return u0 + 0.1 * t;
  }
};

class ExpODE : public MockRHS {
public:
  virtual double dudt(double /* t */, double u) const override { return u; }
  virtual double solution(double t, double u0) const override {
    return u0 * zisa::exp(t);
  }
};

class TSquareODE : public MockRHS {
public:
  virtual double dudt(double t, double /* u */) const override { return t * t; }
  virtual double solution(double t, double u0) const override {
    return u0 + zisa::pow<3>(t) / 3.0;
  }
};

std::shared_ptr<MockRHS> make_mock_rhs(const std::string &desc) {

  if (desc == "constant") {
    return std::make_shared<ConstantRHS>();
  }
  if (desc == "exp") {
    return std::make_shared<ExpODE>();
  }
  if (desc == "t_square") {
    return std::make_shared<TSquareODE>();
  }

  LOG_ERR(string_format("Invalid 'desc'. [%s]", desc.c_str()));
}

class RateOfChangeMock : public RateOfChange {
public:
  RateOfChangeMock(const std::string &desc) : rhs(make_mock_rhs(desc)) {}
  virtual ~RateOfChangeMock() = default;

  /// One update step.
  virtual void
  compute(AllVariables &dudt, const AllVariables &u0, double t) const override {
    for (int_t i = 0; i < u0.size(); ++i) {
      dudt[i] = rhs->dudt(t, u0[i]);
    }
  }

  virtual std::string str() const override { return "RateOfChangeMock"; }

protected:
  std::shared_ptr<MockRHS> rhs;
};

} // namespace zisa

double compute_runge_kutta_error(const std::string &ode_key,
                                 const std::string &solver_key,
                                 double dt) {

  auto dims = zisa::AllVariablesDimensions{30ul, 2ul, 3ul};
  double initial_value = 1.0;

  auto bc = std::make_shared<zisa::DoNothingBC>();
  auto rate_of_change = std::make_shared<zisa::RateOfChangeMock>(ode_key);

  auto solver
      = zisa::make_time_integration(solver_key, rate_of_change, bc, dims);
  auto ode = zisa::make_mock_rhs(ode_key);

  auto u0 = std::make_shared<zisa::AllVariables>(dims);

  for (zisa::int_t i = 0; i < u0->size(); ++i) {
    (*u0)[i] = initial_value;
  }

  double t = 0.0;
  double t_final = 3.0;

  while (t < t_final - 0.5 * dt) {
    u0 = solver->compute_step(u0, t, dt);
    t += dt;
  }

  std::vector<double> l1_err(u0->size());

  for (zisa::int_t i = 0; i < u0->size(); ++i) {
    auto approx = (*u0)[i];
    auto exact = ode->solution(t_final, initial_value);

    l1_err[i] = zisa::abs(approx - exact);
  }

  auto linf_err = std::accumulate(
      l1_err.begin(), l1_err.end(), 0.0, [](double a, double b) {
        return zisa::max(a, b);
      });

  return linf_err;
}

template <class Experiment>
void test_runge_kutta_exact(const std::string &ode_key,
                            const Experiment &experiment) {

  auto &[solver_key, expected_rate, dt_pair] = experiment;
  double dt = 0.01;

  double linf = compute_runge_kutta_error(ode_key, solver_key, dt);

  INFO(string_format(
      "solver = %s, ode = %s", solver_key.c_str(), ode_key.c_str()));
  REQUIRE(linf < 1e-12);
}

template <class Experiment>
void test_runge_kutta_convergence(const std::string &ode_key,
                                  const Experiment &experiment) {

  auto &[solver_key, expected_rate, dt_pair] = experiment;
  auto &[min_rate, max_rate] = expected_rate;
  auto &[dt_coarse, dt_fine] = dt_pair;

  auto linf_coarse = compute_runge_kutta_error(ode_key, solver_key, dt_coarse);
  auto linf_fine = compute_runge_kutta_error(ode_key, solver_key, dt_fine);

  auto rate
      = (log(linf_fine) - log(linf_coarse)) / (log(dt_fine) - log(dt_coarse));

  if (linf_fine > 1e-12) {
    INFO(string_format(
        "solver = %s, ode = %s", solver_key.c_str(), ode_key.c_str()));
    REQUIRE(min_rate <= rate);
    REQUIRE(rate <= max_rate);
  }
}

TEST_CASE("RungeKutta") {

  std::vector<std::string> odes = {"exp", "t_square"};

  using rate_interval_t = std::tuple<double, double>;
  using time_step_pair_t = std::tuple<double, double>;

  std::vector<std::tuple<std::string, rate_interval_t, time_step_pair_t>>
      experiments = {{"ForwardEuler", {0.9, 1.1}, {1e-3, 1e-4}},
                     {"SSP2", {1.9, 2.1}, {1e-3, 1e-4}},
                     {"SSP3", {2.9, 3.1}, {1e-1, 1e-2}},

                     // only second order for non-autonomous ODEs
                     {"Wicker", {1.9, 3.1}, {1e-1, 1e-2}},

                     {"RK4", {3.9, 4.1}, {1e-1, 1e-2}},
                     {"Fehlberg", {4.9, 5.1}, {1e-1, 1e-2}}};

  for (auto &experiment : experiments) {
    test_runge_kutta_exact("constant", experiment);
  }

  for (auto &ode : odes) {
    for (auto &experiment : experiments) {
      test_runge_kutta_convergence(ode, experiment);
    }
  }
}
