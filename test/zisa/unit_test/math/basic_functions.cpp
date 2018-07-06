#include "basic_functions.hpp"
#include <cmath>
#include <cassert>

std::vector<double> convergence_rates(const std::vector<double> &dx,
                                      const std::vector<double> &e) {

  assert(dx.size() == e.size());
  auto r = std::vector<double>(e.size() - 1);

  for (std::size_t i = 0; i < r.size(); ++i) {
    r[i] = (std::log(e[i + 1]) - std::log(e[i]))
           / (std::log(dx[i + 1]) - std::log(dx[i]));
  }

  return r;
}
