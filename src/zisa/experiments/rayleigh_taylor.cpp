#include <zisa/experiments/ic/polytrope_ic.hpp>
#include <zisa/experiments/numerical_experiment_factory.hpp>
#include <zisa/experiments/rayleigh_taylor.hpp>
#include <zisa/parallelization/omp.h>

namespace zisa {

std::shared_ptr<AllVariables> RayleighTaylor::compute_initial_conditions() {
  double amp = params["experiment"]["initial_conditions"]["amplitude"];
  double width = params["experiment"]["initial_conditions"]["width"];
  int n_bumps = params["experiment"]["initial_conditions"]["n_bumps"];

  auto all_variables = compute_initial_conditions(amp, width, n_bumps);
  auto steady_state = compute_initial_conditions(0.0, width, n_bumps);

  auto writer = HDF5SerialWriter(file_name_generator->steady_state_filename);
  save(writer, *steady_state, all_labels<euler_var_t>());

  return all_variables;
}

std::shared_ptr<AllVariables> RayleighTaylor::compute_initial_conditions(
    double amp, double width, int n_bumps) {

  auto dims = choose_all_variable_dims();
  auto all_variables = std::make_shared<AllVariables>(dims);
  const auto &eos = euler->eos;

  const auto &gravity_params = params["euler"]["gravity"];
  auto rhoK = RhoEntropy{gravity_params["rhoC"], gravity_params["K_inner"]};
  auto thetaC = eos.enthalpy_entropy(rhoK);
  auto inner_equilibrium
      = LocalEquilibrium(IsentropicEquilibrium(euler), thetaC, XYZ::zeros());

  double drho = params["experiment"]["initial_conditions"]["drho"];
  double r_crit = gravity_params["r_crit"];
  auto x_ref = XYZ(r_crit * XYZ::unit_vector(0));
  auto rhoE_eq = inner_equilibrium.extrapolate(x_ref);
  auto rhoP = eos.rhoP(rhoE_eq);

  auto theta_inner = eos.enthalpy_entropy(rhoP);
  auto theta_outer = eos.enthalpy_entropy(RhoP{rhoP.rho() + drho, rhoP.p()});

  auto ic_ = PolytropeWithJumpIC(euler, theta_inner, theta_outer, x_ref);
  auto ic = [this, &ic_, amp, width, n_bumps](const auto &x) {
    // Avoid any silent conversion when the return type changes.
    RhoP rhoP = ic_(x);

    double r = zisa::norm(x);
    auto &[rho_eq, p_eq] = rhoP;

    double alpha = zisa::angle(x);

    auto v = amp * exp(-zisa::pow<2>(r / width)) * zisa::sin(n_bumps * alpha);

    double rho = rho_eq;
    double vx = v * zisa::cos(alpha);
    double vy = v * zisa::sin(alpha);
    double p = p_eq;
    double E = euler->energy(rho, vx, vy, p);

    return euler_var_t{rho, rho * vx, rho * vy, 0.0, E};
  };

  auto qr = choose_volume_rule();
  auto &u0 = all_variables->cvars;

  zisa::for_each(triangles(*grid),
                 [&u0, &ic, &qr](int_t i, const Triangle &tri) {
                   u0(i) = average(qr, ic, tri);
                 });

  return all_variables;
}
}
