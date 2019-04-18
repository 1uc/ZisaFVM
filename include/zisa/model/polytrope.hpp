#ifndef ZISA_POLYTROPE_DKOWI_HPP
#define ZISA_POLYTROPE_DKOWI_HPP

#include <zisa/config.hpp>

#include <utility>
#include <zisa/math/cartesian.hpp>
#include <zisa/math/cell_averages.hpp>
#include <zisa/memory/array.hpp>
#include <zisa/model/all_variables.hpp>
#include <zisa/model/gravity.hpp>
#include <zisa/model/radial_poisson_solver.hpp>

namespace zisa {

struct PolytropeParams {
  double rho_center;
  double K_center;
  double G;
  double polytropic_index_n;
  double polytrope_radius;

  double polytropic_gamma() const {
    double n = polytropic_index_n;
    return (n + 1) / n;
  }

  double alpha() const {
    double n = polytropic_index_n;
    return zisa::sqrt((n + 1) * K_center * zisa::pow(rho_center, 1.0 / n - 1.0)
                      / (4 * zisa::pi * G));
  }
};

/// Approximate the Lane-Emden equation on uniform grid.
std::pair<array<double, 1>, array<Cartesian<2>, 1>>
integrate_lane_emden(double r_hat_outer, double n, int_t n_points);

/// Non-uniform version.
std::pair<array<double, 1>, array<Cartesian<2>, 1>>
integrate_lane_emden(array<double, 1> spacing, double r_hat_outer, double n);

void rescale_lane_emden(array<double, 1> &r_hat,
                        array<Cartesian<2>, 1> &x_hat,
                        double rhoC,
                        double alpha,
                        double n);

std::pair<array<double, 1>, array<Cartesian<2>, 1>>
compute_polytrope_profile(const PolytropeParams &params, int_t n_points);

std::pair<array<double, 1>, array<Cartesian<2>, 1>>
compute_polytrope_profile(const PolytropeParams &params,
                          array<double, 1> spacing);

RadialGravity make_general_polytrope_gravity(const PolytropeParams &params,
                                             int_t n_points);

RadialGravity make_general_polytrope_gravity(const PolytropeParams &params,
                                             array<double, 1> spacing);

RadialGravity make_general_polytrope_gravity(double rhoC,
                                             double K,
                                             double G,
                                             double n,
                                             double r_outer,
                                             array<double, 1> spacing);

array<double, 1> compute_polytrope_densities(const Grid &grid,
                                             const PolytropeParams &params);

template <class EOS>
AllVariables make_general_polytrope_profile(const Grid &grid,
                                            const EOS &eos,
                                            const PolytropeParams &params,
                                            int_t quad_deg) {

  auto spacing = zisa::linear_spacing(0.0, 1.0, 2048);
  auto [r, rho_mass] = compute_polytrope_profile(params, spacing);
  auto interpolate = NonUniformLinearInterpolation<Cartesian<2>>(r, rho_mass);

  double K = params.K_center;

  auto cvars = zisa::cell_averages(
      [K, &eos, &interpolate](const XYZ &x) {
        auto rho_mass = interpolate(zisa::norm(x));
        double rho = rho_mass[0];
        auto rhoK = RhoEntropy{rho, K};
        return euler_var_t{rho, 0.0, 0.0, 0.0, eos.internal_energy(rhoK)};
      },
      grid,
      quad_deg);

  return AllVariables{std::move(cvars), GridVariables{}};
}

std::pair<zisa::RadialGravity, std::shared_ptr<zisa::RadialPoissonSolver>>
make_polytrope_self_gravity(const std::shared_ptr<Grid> &grid,
                            const zisa::PolytropeParams &params);

}

#endif // ZISA_POLYTROPE_HPP
