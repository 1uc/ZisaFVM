/* CRTP expressions for Polynomials
 */

#ifndef POLYNOMIAL_EXPR_H_EFZLO
#define POLYNOMIAL_EXPR_H_EFZLO

#include <zisa/math/basic_functions.hpp>
#include <zisa/math/cartesian.hpp>

namespace zisa {

template <class Poly>
class PolynomialCRTP {};

template <class E1, class E2>
class PointwiseSum : public PolynomialCRTP<PointwiseSum<E1, E2>> {
public:
  PointwiseSum(const PolynomialCRTP<E1> &e1, const PolynomialCRTP<E2> &e2)
      : e1(static_cast<const E1 &>(e1)), e2(static_cast<const E2 &>(e2)) {

    static_assert(E1::n_coeffs() == E2::n_coeffs(),
                  "Number of coefficients mismatch.");
    static_assert(E1::n_vars() == E2::n_vars(),
                  "Number of variables mismatch.");
  }

  static constexpr int n_coeffs() { return E1::n_coeffs(); }
  static constexpr int n_vars() { return E1::n_vars(); }

  int degree() const { return zisa::max(e1.degree(), e2.degree()); }
  int n_dims() const { return zisa::max(e1.n_dims(), e2.n_dims()); }

  template <class... Indices>
  double a(Indices &&... indices) const {
    return this->e1.a(std::forward<Indices>(indices)...)
           + this->e2.a(std::forward<Indices>(indices)...);
  }

  template <class... Indices>
  double c(Indices &&... indices) const {
    return (e1.degree() > e2.degree()
                ? e1.c(std::forward<Indices>(indices)...)
                : e2.c(std::forward<Indices>(indices)...));
  }

  const XYZ &x_center() const {
    return (e1.degree() > e2.degree() ? e1.x_center() : e2.x_center());
  }

  double reference_length() const {
    return (e1.degree() > e2.degree() ? e1.reference_length()
                                      : e2.reference_length());
  }

private:
  const E1 &e1;
  const E2 &e2;
};

template <class E1, class E2>
class PointwiseSubtract : public PolynomialCRTP<PointwiseSubtract<E1, E2>> {
public:
  PointwiseSubtract(const PolynomialCRTP<E1> &e1, const PolynomialCRTP<E2> &e2)
      : e1(static_cast<const E1 &>(e1)), e2(static_cast<const E2 &>(e2)) {
    static_assert(E1::n_coeffs() == E2::n_coeffs(),
                  "Number of coefficients mismatch.");
    static_assert(E1::n_vars() == E2::n_vars(),
                  "Number of variables mismatch.");
  }

  static constexpr int n_coeffs() { return E1::n_coeffs(); }
  static constexpr int n_vars() { return E1::n_vars(); }

  int degree() const { return zisa::max(e1.degree(), e2.degree()); }
  int n_dims() const { return zisa::max(e1.n_dims(), e2.n_dims()); }

  template <class... Indices>
  double a(Indices &&... indices) const {
    return e1.a(std::forward<Indices>(indices)...)
           - e2.a(std::forward<Indices>(indices)...);
  }

  template <class... Indices>
  double c(Indices &&... indices) const {
    return (e1.degree() > e2.degree()
                ? e1.c(std::forward<Indices>(indices)...)
                : e2.c(std::forward<Indices>(indices)...));
  }

  const XYZ &x_center() const {
    return (e1.degree() > e2.degree() ? e1.x_center() : e2.x_center());
  }

  double reference_length() const {
    return (e1.degree() > e2.degree() ? e1.reference_length()
                                      : e2.reference_length());
  }

private:
  const E1 &e1;
  const E2 &e2;
};

template <class E>
class PointwiseScale : public PolynomialCRTP<PointwiseScale<E>> {
public:
  PointwiseScale(const PolynomialCRTP<E> &e, double alpha)
      : e(static_cast<const E &>(e)), alpha(alpha) {}

  static constexpr int n_coeffs() { return E::n_coeffs(); }
  static constexpr int n_vars() { return E::n_vars(); }

  int degree() const { return e.degree(); }
  int n_dims() const { return e.n_dims(); }

  template <class... Indices>
  double a(Indices &&... indices) const {
    return alpha * e.a(std::forward<Indices>(indices)...);
  }

  template <class... Indices>
  double c(Indices &&... indices) const {
    return e.c(std::forward<Indices>(indices)...);
  }

  const XYZ &x_center() const { return e.x_center(); }
  double reference_length() const { return e.reference_length(); }

private:
  const E &e;
  double alpha;
};

template <class E1, class E2>
PointwiseSum<E1, E2> operator+(const PolynomialCRTP<E1> &e1,
                               const PolynomialCRTP<E2> &e2) {
  return PointwiseSum<E1, E2>(e1, e2);
}

template <class E1, class E2>
PointwiseSubtract<E1, E2> operator-(const PolynomialCRTP<E1> &e1,
                                    const PolynomialCRTP<E2> &e2) {
  return PointwiseSubtract<E1, E2>(e1, e2);
}

template <class E>
PointwiseScale<E> operator*(double alpha, const PolynomialCRTP<E> &e) {
  return PointwiseScale<E>(e, alpha);
}

template <class E>
PointwiseScale<E> operator/(const PolynomialCRTP<E> &e, double alpha) {
  return PointwiseScale<E>(e, 1.0 / alpha);
}

} // namespace zisa
#endif /* end of include guard */
