/* Shared functionality for ODE rate of change.
 *
 * Authors: Luc Grosheintz <forbugrep@zoho.com>
 *    Date: 2016-09-14
 */
#ifndef RATE_OF_CHANGE_CPP_UYNX1QZJ
#define RATE_OF_CHANGE_CPP_UYNX1QZJ
#include <limits>

#include "zisa/config.hpp"
#include "zisa/math/comparison.hpp"
#include "zisa/model/all_variables.hpp"
#include "zisa/ode/rate_of_change.hpp"
#include "zisa/utils/indent_line.hpp"

namespace zisa {

SumRatesOfChange::SumRatesOfChange(
    std::vector<std::shared_ptr<RateOfChange>> rates_of_change) {
  for (auto &&roc : rates_of_change) {
    add_term(roc);
  }
}

void SumRatesOfChange::compute(AllVariables &tendency,
                               const AllVariables &current_state,
                               double t) {
  for (auto &&roc : rates_of_change) {
    roc->compute(tendency, current_state, t);
  }
}

void SumRatesOfChange::add_term(const std::shared_ptr<RateOfChange> &rate) {
  if (rate != nullptr) {
    rates_of_change.push_back(rate);
  }
}

double
SumRatesOfChange::pick_time_step(const AllVariables &all_variables) const {
  double dt = std::numeric_limits<double>::max();
  return pick_time_step(all_variables, dt);
}

double SumRatesOfChange::pick_time_step(const AllVariables &all_variables,
                                        double dt) const {
  for (auto &&roc : rates_of_change) {
    dt = roc->pick_time_step(all_variables, dt);
  }

  return dt;
}

std::string SumRatesOfChange::str(int indent) const {
  std::stringstream ss;

  for (auto &&roc : rates_of_change) {
    ss << roc->str(indent);
  }

  return ss.str();
}

void SumRatesOfChange::remove_all_terms() { rates_of_change.clear(); }

ZeroRateOfChange::ZeroRateOfChange(double dt_max) : dt_max(dt_max) {}

void ZeroRateOfChange::compute(AllVariables &tendency,
                               const AllVariables &,
                               double) {
  tendency.fill(0);
}

double ZeroRateOfChange::pick_time_step(const AllVariables &) const {
  return dt_max;
}

double ZeroRateOfChange::pick_time_step(const AllVariables &, double dt) const {
  return zisa::min(dt_max, dt);
}

std::string ZeroRateOfChange::str(int indent) const {
  return indent_line(indent, "Zero right-hand side.\n");
}

} // namespace zisa
#endif /* end of include guard: RATE_OF_CHANGE_CPP_UYNX1QZJ */
