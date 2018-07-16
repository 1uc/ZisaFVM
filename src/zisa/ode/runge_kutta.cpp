/* Explicit Runge-Kutta methods.
 */

#include <numeric>

#include <zisa/model/all_variables.hpp>
#include <zisa/ode/runge_kutta.hpp>
#include <zisa/math/comparison.hpp>
#include <zisa/utils/indent_line.hpp>

namespace zisa {
ButcherTableau::ButcherTableau(const std::vector<std::vector<double>> &a,
                               const std::vector<double> &b)
    : n_stages(a.size()) {

  allocate();
  assign(a, b);
}

void ButcherTableau::allocate() {
  this->a.resize(n_stages);
  for (int_t i = 0; i < n_stages; ++i) {
    this->a[i] = array<double, 1>({n_stages}, device);
  }

  this->b = array<double, 1>({n_stages}, device);
  this->c = array<double, 1>({n_stages}, device);
}

void ButcherTableau::assign(const std::vector<std::vector<double>> &a,
                            const std::vector<double> &b) {
  for (int_t i = 0; i < n_stages; ++i) {
    std::copy(a[i].begin(), a[i].end(), this->a[i].begin());
  }

  std::copy(b.begin(), b.end(), this->b.begin());

  for (int_t i = 0; i < n_stages; ++i) {
    c[i] = std::accumulate(a[i].begin(), a[i].end(), 0.0);
  }
}

TendencyBuffers::TendencyBuffers(const AllVariablesDimensions &dims,
                                 int_t n_stages)
    : buffers(n_stages) {
  for (int_t i = 0; i < n_stages; ++i) {
    buffers.at(i) = AllVariables(dims, device);
  }
}

AllVariables &TendencyBuffers::operator[](int_t i) { return buffers.at(i); }

const AllVariables &TendencyBuffers::operator[](int_t i) const {
  return buffers.at(i);
}

std::vector<AllVariables>::iterator TendencyBuffers::begin() {
  return buffers.begin();
}

std::vector<AllVariables>::iterator TendencyBuffers::end() {
  return buffers.end();
}

RungeKutta::RungeKutta(const std::shared_ptr<RateOfChange> &rate_of_change,
                       const std::shared_ptr<BoundaryCondition> &bc,
                       ButcherTableau tableau,
                       const AllVariablesDimensions &dims,
                       double cfl_number)
    : tableau(std::move(tableau)),
      rate_of_change(rate_of_change),
      bc(bc),
      tendency_buffers(dims, this->tableau.n_stages),
      cfl_number(cfl_number) {
  ux = std::make_shared<AllVariables>(dims, device);
}

std::shared_ptr<AllVariables> RungeKutta::compute_step(
    const std::shared_ptr<AllVariables> &u0, double t, double dt) {

  // stage 0
  rate_of_change->compute(tendency_buffers[0], *u0, t);

  // stages 0, ..., s
  for (int_t stage = 1; stage < tableau.n_stages; ++stage) {
    double tx = t + tableau.c[stage] * dt;

    runge_kutta_sum(*ux, *u0, tendency_buffers, tableau.a[stage], dt);
    boundary_condition(*ux, tx);

    rate_of_change->compute(tendency_buffers[stage], *ux, tx);
  }

  runge_kutta_sum(*ux, *u0, tendency_buffers, tableau.b, dt);
  boundary_condition(*ux, t + dt);

  auto tmp = ux;
  ux = u0;
  return tmp;
}

void RungeKutta::boundary_condition(AllVariables &u0, double t) const {
  return bc->apply(u0, t);
}

double RungeKutta::pick_time_step(const AllVariables &all_variables) const {
  return cfl_number * rate_of_change->pick_time_step(all_variables);
}

double RungeKutta::pick_time_step(const AllVariables &all_variables,
                                  double dt) const {
  return zisa::min(dt, pick_time_step(all_variables));
}

std::string RungeKutta::assemble_description(int indent,
                                             const std::string &detail) const {
  return indent_line(indent, detail) + bc->str(indent)
         + rate_of_change->str(indent);
}

void runge_kutta_sum_cpu_impl(AllVariables &u1,
                              const AllVariables &u0,
                              const TendencyBuffers &tendency_buffers,
                              const array<double, 1> &coeffs,
                              double dt) {}

void runge_kutta_sum(AllVariables &u1,
                     const AllVariables &u0,
                     const TendencyBuffers &tendency_buffers,
                     const array<double, 1> &coeffs,
                     double dt) {
  assert(u1.size() == u0.size());

  int_t N = u0.size();
  int_t n_stages = coeffs.shape(0);

  for (int_t i = 0; i < N; ++i) {
    double dudt = 0.0;

    for (int_t stage = 0; stage < n_stages; ++stage) {
      if (coeffs[stage] != 0.0) {
        dudt += coeffs[stage] * tendency_buffers[stage][i];
      }
    }

    u1[i] = u0[i] + dt * dudt;
  }
}

ForwardEuler::ForwardEuler(const std::shared_ptr<RateOfChange> &rate_of_change,
                           const std::shared_ptr<BoundaryCondition> &bc,
                           const AllVariablesDimensions &dims)
    : super(rate_of_change, bc, ButcherTableau({{0.0}}, {1.0}), dims, 0.45) {}

std::string ForwardEuler::str(int indent) const {
  return assemble_description(indent, "Forward Euler (`ForwardEuler`)\n");
}

SSP2::SSP2(const std::shared_ptr<RateOfChange> &rate_of_change,
           const std::shared_ptr<BoundaryCondition> &bc,
           const AllVariablesDimensions &dims)
    : super(rate_of_change,
            bc,
            ButcherTableau({{0.0, 0.0}, {1.0, 0.0}}, {0.5, 0.5}),
            dims,
            0.85) {}

std::string SSP2::str(int indent) const {
  return assemble_description(indent, "SSP 2 (`SSP2`)\n");
}

SSP3::SSP3(const std::shared_ptr<RateOfChange> &rate_of_change,
           const std::shared_ptr<BoundaryCondition> &bc,
           const AllVariablesDimensions &dims)
    : super(
          rate_of_change,
          bc,
          // clang-format off
          ButcherTableau({{0.0, 0.0, 0.0},
                          {1.0, 0.0, 0.0},
                          {0.25, 0.25, 0.0}},
                         {1.0 / 6, 1.0 / 6, 2.0 / 3}),
          // clang-format on
          dims,
          0.85) {}

std::string SSP3::str(int indent) const {
  return assemble_description(indent, "SSP 3 (`SSP3`)\n");
}

Wicker::Wicker(const std::shared_ptr<RateOfChange> &rate_of_change,
               const std::shared_ptr<BoundaryCondition> &bc,
               const AllVariablesDimensions &dims)
    : super(rate_of_change,
            bc,
            // clang-format off
            ButcherTableau(
                {{0.0, 0.0, 0.0},
                 {1.0 / 3, 0.0, 0.0},
                 {0.0, 0.5, 0.0}},
                {0.0, 0.0, 1.0}),
            // clang-format on
            dims,
            0.85) {}

std::string Wicker::str(int indent) const {
  return assemble_description(indent, "Wicker (`Wicker`)\n");
}

RK4::RK4(const std::shared_ptr<RateOfChange> &rate_of_change,
         const std::shared_ptr<BoundaryCondition> &bc,
         const AllVariablesDimensions &dims)
    : super(rate_of_change,
            bc,
            // clang-format off
                ButcherTableau({{0.0, 0.0, 0.0, 0.0},
                                {0.5, 0.0, 0.0, 0.0},
                                {0.0, 0.5, 0.0, 0.0},
                                {0.0, 0.0, 1.0, 0.0}},
                                {1.0 / 6, 1.0 / 3, 1.0 / 3, 1.0 / 6}),
            // clang-format on
            dims,
            0.85) {}

std::string RK4::str(int indent) const {
  return assemble_description(indent, "The Runge Kutta (`RK4`)\n");
}

Fehlberg::Fehlberg(const std::shared_ptr<RateOfChange> &rate_of_change,
                   const std::shared_ptr<BoundaryCondition> &bc,
                   const AllVariablesDimensions &dims)
    : super(rate_of_change,
            bc,
            // clang-format off
            ButcherTableau(
                {{0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
                 {0.25, 0.0, 0.0, 0.0, 0.0, 0.0},
                 {3.0 / 32.0, 9.0 / 32.0, 0.0, 0.0, 0.0, 0.0},
                 {1932.0 / 2197.0,
                  -7200.0 / 2197.0,
                  7296.0 / 2197.0,
                  0.0,
                  0.0,
                  0.0},
                 {439.0 / 216.0, -8.0, 3680.0 / 513.0, -845.0 / 4104, 0.0, 0.0},
                 {-8.0 / 27.0,
                  2.0,
                  -3544.0 / 2565.0,
                  1859.0 / 4104.0,
                  -11.0 / 40.0,
                  0.0}},

                {16.0 / 135.0,
                 0.0,
                 6656.0 / 12825.0,
                 28561.0 / 56430.0,
                 -9.0 / 50.0,
                 2.0 / 55.0}),
            // clang-format on
            dims,
            0.85) { /* otherwise empty constructor */
}

std::string Fehlberg::str(int indent) const {
  return assemble_description(indent,
                              "Fehlberg RK method, sixth order part.\n");
}

} // namespace zisa