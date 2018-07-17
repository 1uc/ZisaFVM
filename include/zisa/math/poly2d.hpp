/*
 *
 */

#ifndef POLY2D_H_HRE31
#define POLY2D_H_HRE31

#include <zisa/config.hpp>
#include <zisa/math/cartesian.hpp>
#include <zisa/math/polynomial_expr.hpp>
#include <zisa/io/format_as_list.hpp>

namespace zisa {

constexpr ANY_DEVICE_INLINE int_t poly_dof(int deg) {
  return int_t((deg + 1) * (deg + 2) / 2);
}

/// Highest degree of a polynomial with '<=n_coeffs' coefficients.
int poly_degree(int_t n_coeffs);

/// Compute the linear index of the coefficient of a polynomial.
/** Example:
 *    p(x, y) = \sum_{k,l} a(k,l) * x^k*y^l
 *
 *    (0, 0) -> 0,
 *    (1, 0) -> 1, (0, 1) -> 2,
 *    (2, 0) -> 3, (1, 1) -> 4, (0, 2) -> 5,
 *    ...
 *
 *  Input:
 *    @param k  -- slow index
 *    @param l  -- fast index, i.e. for fixed (k+l) consecutive values of l
 *                 result in consecutive values of the linear index.
 **/
constexpr ANY_DEVICE_INLINE int_t poly_index(int a, int b) {
  auto n = a + b;
  return int_t(((n + 1) * n) / 2 + b);
}

template <int MAX_DEGREE>
class Poly2D : public PolynomialCRTP<Poly2D<MAX_DEGREE>> {
public:
  static constexpr int max_degree() { return MAX_DEGREE; }
  int degree() const { return degree_; }

  static constexpr int_t idx(int k, int l) { return poly_index(k, l); }
  static constexpr int_t dof(int deg) { return poly_dof(deg); }
  static constexpr int_t n_coeffs() { return dof(max_degree()); }

public:
  Poly2D(void) : degree_(0) { std::fill(coeffs, coeffs + n_coeffs(), 0.0); }

  Poly2D(const std::initializer_list<double> &coeffs_list,
         const std::initializer_list<double> &moments_list)
      : degree_(poly_degree(coeffs_list.size())) {

    assert(coeffs_list.size() <= n_coeffs());
    assert(coeffs_list.size() == moments_list.size());

    // no excess arguments.
    assert(poly_dof(degree()) == coeffs_list.size());

    std::copy(coeffs_list.begin(), coeffs_list.end(), coeffs);
    std::copy(moments_list.begin(), moments_list.end(), moments);
  }

  template <class E>
  Poly2D(const PolynomialCRTP<E> &e_) {
    deep_copy(static_cast<const E &>(e_));
  }

  template <int D>
  void operator=(const Poly2D<D> &other) {
    deep_copy(other);
  }

  template <class E>
  void operator=(const PolynomialCRTP<E> &e_) {
    const E &e = static_cast<const E &>(e_);

    deep_copy(e);
  }

  template <class E>
  void deep_copy(const PolynomialCRTP<E> &e_) {
    const E &e = static_cast<const E &>(e_);

    static_assert(E::max_degree() <= max_degree(),
                  "Polynomial degree mismatch.");

    degree_ = e.degree();

    int deg = degree();
    for (int d = 0; d <= deg; ++d) {
      for (int k = 0; k <= d; ++k) {
        int l = d - k;
        this->a(k, l) = e.a(k, l);
        this->c(k, l) = e.c(k, l);
      }
    }
  }

  double operator()(const XY &xy) const {
    auto x = xy[0];
    auto y = xy[1];

    auto px = a(0, 0);

    for (int d = 1; d <= degree(); ++d) {
      for (int k = 0; k <= d; ++k) {
        int l = d - k;
        px += a(k, l) * (std::pow(x, k) * std::pow(y, l) - c(k, l));
      }
    }

    return px;
  }

  double &a(int i, int j) { return coeffs[idx(i, j)]; }
  double &c(int i, int j) { return moments[idx(i, j)]; }

  double a(int i, int j) const {
    assert(i >= 0);
    assert(j >= 0);

    if (i + j <= degree()) {
      return coeffs[idx(i, j)];
    } else {
      return 0.0;
    }
  }

  double c(int i, int j) const {
    assert(i >= 0);
    assert(j >= 0);

    LOG_ERR_IF(i + j > degree(), "Don't have the moments for this.");

    return moments[idx(i, j)];
  }

protected:
  int degree_;
  double coeffs[n_coeffs()];
  double moments[n_coeffs()];

  template <int DEG>
  friend std::ostream &operator<<(std::ostream &os,
                                  const Poly2D<DEG> &poly2d);
};

template <int MAX_DEGREE>
std::ostream &operator<<(std::ostream &os, const Poly2D<MAX_DEGREE> &poly2d) {
  os << format_as_list(poly2d.coeffs, poly_dof(poly2d.degree())) << "\n";
  os << format_as_list(poly2d.moments, poly_dof(poly2d.degree())) << "\n";

  return os;
}

} // namespace zisa
#endif /* end of include guard */
