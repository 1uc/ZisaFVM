#include <zisa/experiments/ic/polytrope_ic.hpp>
#include <zisa/experiments/rayleigh_taylor.hpp>

namespace zisa {

std::shared_ptr<AllVariables> RayleighTaylor::compute_initial_conditions() {
  double amp = params["experiment"]["initial_conditions"]["amplitude"];
  double amp_noise = params["experiment"]["initial_conditions"]["A_noise"];
  double width = params["experiment"]["initial_conditions"]["width"];
  int n_bumps = params["experiment"]["initial_conditions"]["n_bumps"];

  auto all_variables
      = compute_initial_conditions(amp, amp_noise, width, n_bumps);
  auto steady_state = compute_initial_conditions(0.0, 0.0, width, n_bumps);

  auto vis = choose_visualization();
  vis->steady_state(*steady_state);

  return all_variables;
}

int_t RayleighTaylor::choose_n_avars() { return 1; }

std::shared_ptr<AllVariables> RayleighTaylor::compute_initial_conditions(
    double amp, double amp_noise, double width, int n_bumps) {

  auto dims = choose_all_variable_dims();
  auto all_variables = std::make_shared<AllVariables>(dims);
  auto grid = choose_grid();
  const auto &local_eos = choose_local_eos();
  const auto &eos = (*local_eos)(0);

  const auto &gravity_params = params["euler"]["gravity"];
  auto rhoK = RhoEntropy{gravity_params["rhoC"], gravity_params["K"]};
  auto thetaC = eos->enthalpy_entropy(rhoK);
  auto inner_equilibrium = LocalEquilibrium(
      IsentropicEquilibrium(eos, gravity), thetaC, XYZ::zeros());

  const auto &ic_params = params["experiment"]["initial_conditions"];
  double amp_shape = ic_params["A_shape"];
  double drho = ic_params["drho"];
  double r_crit = ic_params["r_crit"];
  auto x_ref = XYZ(r_crit * XYZ::unit_vector(0));
  auto rhoE_eq = inner_equilibrium.extrapolate(x_ref);
  auto rhoP = eos->rhoP(rhoE_eq);

  double r_noise = ic_params["r_noise"];

  auto theta_inner = eos->enthalpy_entropy(rhoP);
  auto theta_outer = eos->enthalpy_entropy(RhoP{rhoP.rho() + drho, rhoP.p()});

  auto ic_ = PolytropeWithJumpIC(
      eos, gravity, theta_inner, theta_outer, x_ref, amp_shape, n_bumps);
  auto ic = [this, &eos, &ic_, amp, width, n_bumps, r_crit, r_noise, amp_noise](
                const auto &x) {
    // Avoid any silent conversion when the return type changes.
    RhoP rhoP = ic_(x);

    double r = zisa::norm(x);
    auto &[rho_eq, p_eq] = rhoP;

    double alpha = zisa::angle(x);

    auto v = amp * exp(-zisa::pow<2>((r - r_crit) / width))
             * zisa::sin(n_bumps * alpha);

    double noise_width = 2.5e-02;

    double rho = rho_eq;
    double vx = v * zisa::cos(alpha);
    double vy = v * zisa::sin(alpha);
    double p = p_eq;

    for (int k = 0; k < n_bumps; ++k) {
      double beta = k * 2.0 * pi / n_bumps;
      auto x_noise
          = XYZ{r_noise * zisa::cos(beta), r_noise * zisa::sin(beta), 0.0};

      auto drho = amp_noise
                  * exp(-zisa::pow<2>(zisa::norm(x - x_noise) / noise_width))
                  * zisa::sin2pi(0.25 * (r - r_noise) / noise_width);

      rho += drho;
    }

    double E = total_energy(*eos, rho, vx, vy, p);
    return euler_var_t{rho, rho * vx, rho * vy, 0.0, E};
  };

  auto tracer_ic = [&ic_, &ic](const auto &x) {
    return ic(x)[0] * (ic_.is_inner(x) ? 1.0 : 10.0);
  };

  auto qr = choose_volume_rule();
  auto &u0 = all_variables->cvars;
  auto &q0 = all_variables->avars;

  zisa::for_each(
      serial_policy{},
      triangles(*grid),
      [&q0, &u0, &ic, &tracer_ic, &qr](int_t i, const Triangle &tri) {
        u0(i) = average(qr, ic, tri);
        q0(i, 0) = average(qr, tracer_ic, tri);
      });

  return all_variables;
}

std::function<bool(const Grid &grid, int_t i)>
RayleighTaylor::boundary_mask() const {

  return [](const Grid &grid, int_t i) {
    // FIXME hard-coded value.
    double r_crit = 0.6;
    return zisa::norm(grid.cell_centers[i]) > r_crit;
  };
}

}
