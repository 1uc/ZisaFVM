#include <zisa/math/denormalized_rule.hpp>

#include <numeric>

namespace zisa {
DenormalizedRule::DenormalizedRule(int_t n_points)
    : weights(n_points), points(n_points) {

  volume = std::accumulate(weights.begin(), weights.end(), 0.0);
}

bool operator!=(const DenormalizedRule &a, const DenormalizedRule &b) {
  return !(a == b);
}

bool operator==(const DenormalizedRule &a, const DenormalizedRule &b) {
  return (a.weights == b.weights) && (a.points == b.points)
         && (a.volume == b.volume);
}
}