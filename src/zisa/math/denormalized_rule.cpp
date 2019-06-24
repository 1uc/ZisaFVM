#include <zisa/math/denormalized_rule.hpp>

namespace zisa {
DenormalizedRule::DenormalizedRule(int_t n_points)
    : weights(n_points), points(n_points) {}
}