/*
 *
 * Authors: Luc Grosheintz <forbugrep@zoho.com>
 *    Date: 2017-07-16
 */

#ifndef GAUSS_LEGENDRE_H_S48VB
#define GAUSS_LEGENDRE_H_S48VB

#include <type_traits>

#include <zisa/math/fourier_series.hpp>
#include <zisa/math/newton.hpp>

namespace zisa {
// Forward declarations
template <int n>
struct GaussLegendre;

template <int n>
void fourier_newton(GaussLegendre<n> &gl);

/// Normalized to [-1, 1].
template <int n>
struct GaussLegendre {
  GaussLegendre() { fourier_newton(*this); }

  double points[n];
  double weights[n];

  static inline constexpr int n_qp() { return n; };
};

/// Find Gauss-Legendre quadrature point by root finding.
/** The quadrature points of Gauss-Legendre are the roots of the
 * Legendre polynomials.
 *
 * This algorithm utilized the fact the roots of the Fourier
 * series of the Legendre-polynomial are nearly evenly spaced. Hence,
 * a good initial guess for Newton's method exists and the roots are
 * found using Newtons method.
 *
 * References:
 *   [0]: Belousov, Legendre Polynomials: Mathematical Tables Series, 1962.
 *   [1]: Swarztrauber, SIAM J. Sci. Comp. 2002
 *   [2]: Hale, Townsend, SIAM J. Sci. Comp. 2013
 */
class FourierNewton {
public:
  template <int n>
  void operator()(GaussLegendre<n> &gl) {
    auto p = fourier_poly<n>();
    auto dp = p.derivative();

    quadrature_points(gl, p, dp);
    quadrature_weights(gl, p, dp);
    convert_theta_to_x(gl);
  }

  template <int n>
  void quadrature_points(GaussLegendre<n> &gl,
                         const FourierSeries<n> &p,
                         const FourierSeriesDerivative<n> &dp) {
    bool is_even = n % 2 == 0;

    if (is_even) {
      quadrature_points_even(gl, p, dp);
    } else {
      quadrature_points_odd(gl, p, dp);
    }
  }

  template <int n>
  void quadrature_points_odd(GaussLegendre<n> &gl,
                             const FourierSeries<n> &p,
                             const FourierSeriesDerivative<n> &dp) {
    gl.points[n / 2] = zisa::pi / 2;

    for (int k = 1; k <= n / 2; ++k) {
      double theta_guess
          = (k < 2 ? zisa::pi / 2 - k * zisa::pi / n
                   : 2.0 * gl.points[n / 2 + k - 1] - gl.points[n / 2 + k - 2]);

      gl.points[n / 2 + k] = newton(p, dp, theta_guess, 1e-12);
    }
  }

  template <int n>
  void quadrature_points_even(GaussLegendre<n> &gl,
                              const FourierSeries<n> &p,
                              const FourierSeriesDerivative<n> &dp) {
    for (int k = 0; k < n / 2; ++k) {
      double theta_guess
          = (k < 2 ? zisa::pi / 2 - (k + 0.5) * zisa::pi / n
                   : 2.0 * gl.points[n / 2 + k - 1] - gl.points[n / 2 + k - 2]);

      gl.points[n / 2 + k] = newton(p, dp, theta_guess, 1e-12);
    }
  }

  template <int n>
  void quadrature_weights(GaussLegendre<n> &gl,
                          const FourierSeries<n> &,
                          const FourierSeriesDerivative<n> &dp) {
    bool is_even = n % 2 == 0;
    bool is_odd = !is_even;

    if (is_odd) {
      gl.weights[n / 2] = weight(dp, gl.points[n / 2]);
    }

    for (int k = 0; k < n / 2; ++k) {
      int i_upper = (is_even ? n / 2 + k : n / 2 + k + 1);
      int i_lower = n / 2 - 1 - k;

      gl.weights[i_upper] = weight(dp, gl.points[i_upper]);
      gl.weights[i_lower] = gl.weights[i_upper];
    }
  }

  template <int n>
  double weight(const FourierSeriesDerivative<n> &dp, double theta) {
    return (2 * n + 1) / zisa::pow<2>(dp(theta));
  }

  template <int n>
  void convert_theta_to_x(GaussLegendre<n> &gl) {
    bool is_even = n % 2 == 0;
    bool is_odd = !is_even;

    if (is_odd) {
      gl.points[n / 2] = zisa::cos(gl.points[n / 2]);
    }

    for (int k = 0; k < n / 2; ++k) {
      int i_upper = (is_even ? n / 2 + k : n / 2 + k + 1);
      int i_lower = n / 2 - 1 - k;

      gl.points[i_upper] = zisa::cos(gl.points[i_upper]);
      gl.points[i_lower] = -gl.points[i_upper];
    }
  }

  template <int n>
  FourierSeries<n> fourier_poly() {
    constexpr bool is_even = n % 2 == 0;

    FourierSeries<n> p;
    for (int i = 0; i <= n; ++i) {
      p.at(i) = 0.0;
    }

    p.at(n) = a_nn(n);
    for (int l = 2; l <= n; l += 2) {
      double factor
          = double((l - 1) * (2 * n - l + 2)) / double(l * (2 * n - l + 1));
      p.at(n - l) = factor * p.at(n - l + 2);
    }

    if (is_even) {
      p.at(0) *= 0.5;
    }

    return p;
  }

protected:
  double a_nn(int n) {
    double a = zisa::sqrt(2.0);

    for (int i = 1; i <= n; ++i) {
      a = zisa::sqrt(1.0 - 1.0 / (4.0 * zisa::pow<2>(i))) * a;
    }

    return a;
  }
};

template <int n>
void fourier_newton(GaussLegendre<n> &gl) {
  auto fn = FourierNewton();
  fn(gl);
}

template <int n_qp>
const GaussLegendre<n_qp> &make_gauss_legendre() {
  static GaussLegendre<n_qp> the_one;
  return the_one;
}

} // namespace zisa
#endif /* end of include guard */
