#ifndef ZISA_LINEAR_INTERPOLATION_EKOSI_HPP
#define ZISA_LINEAR_INTERPOLATION_EKOSI_HPP
#include <zisa/config.hpp>

#include <algorithm>
#include <memory>
#include <zisa/memory/array.hpp>

namespace zisa {
template <class T>
class NonUniformLinearInterpolation {
public:
  NonUniformLinearInterpolation() = default;
  NonUniformLinearInterpolation(array<double, 1> points, array<T, 1> values)
      : points(std::move(points)), values(std::move(values)) {}

  T operator()(double r) const {
    int_t i = index(r);
    double alpha = (r - points(i)) / (points(i + 1) - points(i));
    return T((1 - alpha) * values[i] + alpha * values[i + 1]);
  }

  T derivative(double r) const {
    int_t i = index(r);
    return T((values[i + 1] - values[i]) / (points(i + 1) - points(i)));
  }

  int_t index(double r) const {
    auto first_past = std::lower_bound(points.cbegin(), points.cend(), r);

    auto i = int_t(first_past - points.begin());
    if (i == 0) {
      return i;
    }
    assert(i < values.size());
    return i - 1;
  }

public:
  array<double, 1> points;
  array<T, 1> values;
};

}

#endif // ZISA_LINEAR_INTERPOLATION_HPP
