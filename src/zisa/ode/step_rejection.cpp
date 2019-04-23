#include <zisa/ode/step_rejection.hpp>

namespace zisa {
bool RejectNothing::is_good() const { return true; }

double RejectNothing::pick_time_step(double dt) const { return dt; }

void RejectNothing::check(const AllVariables &, const AllVariables &) const {
  // do nothing
}

bool RejectLargeDensityChange::is_good() const { return is_good_; }

void RejectLargeDensityChange::check(const AllVariables &all_vars0,
                                     const AllVariables &all_vars1) const {
  const auto &u0 = all_vars0.cvars;
  const auto &u1 = all_vars1.cvars;

  int_t n_cells = u0.shape(0);

  bool is_this_step_good = true;
#pragma omp parallel for reduction(&& : is_this_step_good)
  for (int_t i = 0; i < n_cells; ++i) {
    if (zisa::abs(u0(i, 0) - u1(i, 0)) > drho_crit_rel * u0(i, 0)) {
      is_this_step_good = false;
    }
  }

  if (!is_this_step_good) {
    LOG_WARN("Shrinking step size.");
    is_good_ = false;
    current_factor /= growth_factor;
  } else {
    if (is_good_) {
      // 2 consecutive steps succeeded.
      current_factor *= growth_factor;
    }
    is_good_ = true;
  }

  current_factor = zisa::max(minimum_factor, zisa::min(current_factor, 1.0));
}

double RejectLargeDensityChange::pick_time_step(double dt) const {
  return current_factor * dt;
}

RejectLargeDensityChange::RejectLargeDensityChange(double drho_crit_rel,
                                                   int max_exponent)
    : drho_crit_rel(drho_crit_rel),
      minimum_factor(zisa::pow(growth_factor, -max_exponent)) {}
}