#ifndef RATE_OF_CHANGE_CPP_UYNX1QZJ
#define RATE_OF_CHANGE_CPP_UYNX1QZJ

#include <zisa/config.hpp>
#include <zisa/model/all_variables.hpp>
#include <zisa/ode/rate_of_change.hpp>
#include <zisa/utils/timer.hpp>

namespace zisa {

SumRatesOfChange::SumRatesOfChange(
    const std::vector<std::shared_ptr<RateOfChange>> &rates_of_change) {
  for (auto &&roc : rates_of_change) {
    add_term(roc);
  }
}

void SumRatesOfChange::compute(AllVariables &tendency,
                               const AllVariables &current_state,
                               double t) const {
  for (auto &&roc : rates_of_change) {
    roc->compute(tendency, current_state, t);
  }
}

void SumRatesOfChange::add_term(const std::shared_ptr<RateOfChange> &rate) {
  if (rate != nullptr) {
    rates_of_change.push_back(rate);
  }
}

std::string SumRatesOfChange::str() const {
  std::stringstream ss;
  bool is_first = true;

  for (auto &&roc : rates_of_change) {
    ss << (is_first ? "" : "\n") << roc->str();
    is_first = false;
  }

  return ss.str();
}

void SumRatesOfChange::remove_all_terms() { rates_of_change.clear(); }

void ZeroRateOfChange::compute(AllVariables &tendency,
                               const AllVariables & /* current_state */,
                               double /* t */) const {
  zisa::fill(tendency.cvars, 0.0);
  zisa::fill(tendency.avars, 0.0);
}

std::string ZeroRateOfChange::str() const { return "Zero right-hand side."; }

} // namespace zisa
#endif /* end of include guard: RATE_OF_CHANGE_CPP_UYNX1QZJ */
