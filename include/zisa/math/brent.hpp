#ifndef ZISA_BRENT_XEDSW_HPP
#define ZISA_BRENT_XEDSW_HPP

#include <zisa/config.hpp>

#include <array>
#include <utility>
#include <zisa/math/basic_functions.hpp>

namespace zisa {

inline double
inverse_quadratic_interpolation_step(const std::array<double, 2> &a_fa,
                                     const std::array<double, 2> &b_fb,
                                     const std::array<double, 2> &c_fc) {

  auto inv = [](const std::array<double, 2> &a_fa,
                const std::array<double, 2> &b_fb,
                const std::array<double, 2> &c_fc) {
    auto [a, fa] = a_fa;
    auto [b, fb] = b_fb;
    auto [c, fc] = c_fc;

    return a * fb * fc / ((fa - fb) * (fa - fc));
  };

  return inv(a_fa, b_fb, c_fc) + inv(b_fb, a_fa, c_fc) + inv(c_fc, a_fa, b_fb);
}

inline double secant_method_step(const std::array<double, 2> &a_fa,
                                 const std::array<double, 2> &b_fb) {

  auto [a, fa] = a_fa;
  auto [b, fb] = b_fb;

  return b - fb * (b - a) / (fb - fa);
}

template <class F>
std::pair<double, double> find_bracket(const F &f,
                                       const std::array<double, 2> &range_guess,
                                       const std::array<double, 2> &hard_bounds,
                                       int max_iter) {
  auto [x_min, x_max] = hard_bounds;
  auto [a, b] = range_guess;

  if (a > b) {
    zisa::swap(a, b);
  }

  auto shrink = [x_min = x_min, a0 = a, b0 = b](double a, int i) {
    if (zisa::isreal(x_min)) {
      return zisa::max(a - pow(b0 - a0, i), 0.5 * (a + x_min));
    } else {
      return a - pow(b0 - a0, i);
    }
  };

  auto grow = [x_max = x_max, a0 = a, b0 = b](double b, int i) {
    if (zisa::isreal(x_max)) {
      return zisa::min(b + pow(b0 - a0, i), 0.5 * (b + x_max));
    } else {
      return b + pow(b0 - a0, i);
    }
  };

  double fa = f(a);
  double fb = f(b);

  if (fa * fb < 0.0) {
    return {a, b};
  }

  for (int i = 0; i < max_iter; ++i) {
    a = shrink(a, i);
    fa = f(a);

    if (fa * fb < 0.0) {
      break;
    }

    b = grow(b, i);
    fb = f(b);

    if (fa * fb < 0.0) {
      break;
    }
  }

  return {a, b};
}

/// Brent's algorithm for solving non-linear equations.
/**
 *  References:
 *   [1]: https://en.wikipedia.org/wiki/Brent's_method
 *
 *  Parameters:
 *       f : Callable with signature `double f(double)`.
 *
 *    a, b : Non-empty interval in which the root must lie, i.e.
 *              x_star \in [a, b] or x_star \in [b, a].
 */

template <class F>
double brent(const F &f, double a, double b, double atol, int max_iter) {
  double delta = 1e-6;
  auto is_outside_of = [](double x, double a, double b) {
    return (a < b) ? (x < a || x > b) : (x < b || x > a);
  };

  double fa = f(a);
  double fb = f(b);

  LOG_ERR_IF(fa * fb >= 0,
             string_format(
                 "This is not a bracket: (%e, %e), (%e, %e).", a, fa, b, fb));

  if (zisa::abs(fa) < zisa::abs(fb)) {
    zisa::swap(a, b);
    zisa::swap(fa, fb);
  }

  double d = a;
  double c = a;
  double s = a;

  double fc = f(c);
  double fs = f(s);

  bool used_bisection = true;

  for (int i = 0; i < max_iter; ++i) {
    if (almost_equal(fb, 0.0, atol)) {
      return b;
    }

    if (almost_equal(fs, 0.0, atol)) {
      return s;
    }

    if (almost_equal(a, b, atol)) {
      return a;
    }

    if (!almost_equal(fa, fc, atol) || !almost_equal(fb, fc, atol)) {
      s = inverse_quadratic_interpolation_step({a, fa}, {b, fb}, {c, fc});
    } else {
      s = secant_method_step({a, fa}, {b, fb});
    }

    if (is_outside_of(s, 0.25 * (3 * a + b), b)
        || (used_bisection && zisa::abs(s - b) >= 0.5 * zisa::abs(b - c))
        || (!used_bisection && zisa::abs(s - b) >= 0.5 * zisa::abs(c - d))
        || (used_bisection && zisa::abs(b - c) < delta)
        || (!used_bisection && zisa::abs(c - d) < delta)) {

      s = 0.5 * (a + b);
      used_bisection = true;
    } else {
      used_bisection = false;
    }

    fs = f(s);
    d = c;
    c = b;
    fc = fb;

    if (fa * fs < 0.0) {
      b = s;
      fb = f(b);
    } else {
      a = s;
      fa = f(a);
    }

    if (zisa::abs(fa) < zisa::abs(fb)) {
      zisa::swap(a, b);
      zisa::swap(fa, fb);
    }
  }

  LOG_ERR("Failed to converge.");
}

}
#endif // ZISA_BRENT_HPP
