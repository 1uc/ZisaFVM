#include <zisa/math/linear_spacing.hpp>
#include <zisa/model/polytrope.hpp>
#include <zisa/ode/runge_kutta.hpp>

namespace zisa {
std::pair<array<double, 1>, array<Cartesian<2>, 1>>
integrate_lane_emden(double r_hat_outer, double n, int_t n_points) {
  auto spacing = zisa::linear_spacing(0, 1.0, n_points);
  return integrate_lane_emden(spacing, r_hat_outer, n);
}

std::pair<array<double, 1>, array<Cartesian<2>, 1>>
integrate_lane_emden(array<double, 1> spacing, double r_hat_outer, double n) {
  // Numerically solve:
  //
  // d \rho_hat / d \r_hat = - \m_hat / \r_hat^2
  // d \m_hat / d \r_hat = \r_hat^2 \rho_hat^n
  //
  // dx / dt = f(t, x)

  using X = Cartesian<2>;

  auto f = [n](double t, const X &x) {
    double eps = std::numeric_limits<double>::min();
    return X{-x[1] / (pow<2>(t) + eps), pow<2>(t) * zisa::pow(x[0], n)};
  };

  auto n_points = spacing.size();
  auto x = array<X, 1>(n_points);
  auto rk = StaticRungeKutta<X>(f, "rk4");
  x[0] = X{1.0, 0.0};

  spacing[0] = 0.0;

  for (int_t i = 0; i < n_points - 1; ++i) {
    spacing[i + 1] = r_hat_outer * spacing[i + 1];
    double dr_hat = spacing[i + 1] - spacing[i];

    x[i + 1] = rk(x[i], spacing[i], dr_hat);
  }
  spacing[n_points - 1] = r_hat_outer;

  return std::pair<array<double, 1>, array<Cartesian<2>, 1>>(std::move(spacing),
                                                             std::move(x));
}

void rescale_lane_emden(array<double, 1> &r_hat,
                        array<Cartesian<2>, 1> &x_hat,
                        double rhoC,
                        double alpha,
                        double n) {

  for (auto &x : x_hat) {
    const auto &[rho_hat, m_hat] = x;

    x[0] = rhoC * zisa::pow(rho_hat, n);
    x[1] = 4 * pi * zisa::pow<3>(alpha) * rhoC * m_hat;
  }

  for (auto &r : r_hat) {
    r = alpha * r;
  }
}

std::pair<array<double, 1>, array<Cartesian<2>, 1>>
compute_polytrope_profile(double r_outer,
                          double rhoC,
                          double alpha,
                          double polytropic_index,
                          array<double, 1> spacing) {
  double r_hat_outer = r_outer / alpha;

  auto [r_hat, x_hat]
      = integrate_lane_emden(std::move(spacing), r_hat_outer, polytropic_index);
  rescale_lane_emden(r_hat, x_hat, rhoC, alpha, polytropic_index);

  return std::pair<array<double, 1>, array<Cartesian<2>, 1>>(std::move(r_hat),
                                                             std::move(x_hat));
}

std::pair<array<double, 1>, array<Cartesian<2>, 1>>
compute_polytrope_profile(const PolytropeParams &params, int_t n_points) {

  auto spacing = zisa::linear_spacing(0.0, 1.0, n_points);
  return compute_polytrope_profile(params, std::move(spacing));
}

std::pair<array<double, 1>, array<Cartesian<2>, 1>>
compute_polytrope_profile(const PolytropeParams &params,
                          array<double, 1> spacing) {
  double r_outer = params.polytrope_radius;
  double rhoC = params.rho_center;
  double n = params.polytropic_index_n;
  double alpha = params.alpha();

  return compute_polytrope_profile(r_outer, rhoC, alpha, n, std::move(spacing));
}

RadialGravity make_general_polytrope_gravity(const PolytropeParams &params,
                                             int_t n_points) {
  auto spacing = zisa::linear_spacing(0.0, 1.0, n_points);
  return make_general_polytrope_gravity(params, std::move(spacing));
}

RadialGravity make_general_polytrope_gravity(const PolytropeParams &params,
                                             array<double, 1> spacing) {

  double G = params.G;
  double rhoC = params.rho_center;
  double alpha = params.alpha();
  double n = params.polytropic_index_n;
  double r_outer = params.polytrope_radius;
  int_t n_points = spacing.size();

  auto [r, x]
      = compute_polytrope_profile(r_outer, rhoC, alpha, n, std::move(spacing));

  auto phi = array<double, 1>(n_points);
  phi[0] = 0.0;
  for (int_t i = 1; i < n_points; ++i) {
    double m = x(i)[1];
    double r2 = zisa::pow<2>(r(i));
    double dr = r(i) - r(i - 1);
    phi[i] = phi[i - 1] + G * m / r2 * dr;
  }

  return RadialGravity(std::move(r), std::move(phi));
}

std::pair<zisa::RadialGravity, std::shared_ptr<zisa::RadialPoissonSolver>>
make_polytrope_self_gravity(const std::shared_ptr<Grid> &grid,
                            const zisa::PolytropeParams &params) {
  auto G = params.G;
  auto [gravity, poisson_solver] = zisa::make_radial_poisson_solver(grid, G);

  // set gravity
  {
    auto eos = zisa::IdealGasEOS(params.polytropic_gamma(), -1.0);
    auto all_variables
        = zisa::make_general_polytrope_profile(*grid, eos, params);
    poisson_solver->update(gravity, all_variables);
  }

  return {gravity, poisson_solver};
}

}