/* Two dimensional polynomial used for WENO reconstruction.
 */
#ifndef POLY2D_DECL_H_83W71
#define POLY2D_DECL_H_83W71

#include <zisa/config.hpp>
#include <zisa/io/format_as_list.hpp>
#include <zisa/math/cartesian.hpp>
#include <zisa/math/polynomial_expr.hpp>
#include <zisa/memory/array.hpp>

namespace zisa {

/// Number of degree of freedom in a `deg` degree 2D polynomial.
constexpr ANY_DEVICE_INLINE int_t poly_dof(int deg);

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
constexpr ANY_DEVICE_INLINE int_t poly_index(int k, int l);

template <int MAX_DEGREE, int NVARS>
class Poly2D : public PolynomialCRTP<Poly2D<MAX_DEGREE, NVARS>> {
public:
  static constexpr int max_degree();
  int degree() const;

  static constexpr int_t idx(int k, int l);
  static constexpr int_t dof(int deg);
  static constexpr int_t n_coeffs();
  static constexpr int_t n_vars();

public:
  Poly2D();

  Poly2D(int degree,
         const array<double, 1> &moments,
         const XY &x_center,
         double reference_length);

  Poly2D(int degree,
         const std::initializer_list<double> &moments_list,
         const XY &x_center,
         double reference_length);

  Poly2D(const std::initializer_list<double> &coeffs_list,
         const std::initializer_list<double> &moments_list,
         const XY &x_center,
         double reference_length);

  template <class E>
  Poly2D(const PolynomialCRTP<E> &e_);

  template <class E>
  void operator=(const PolynomialCRTP<E> &e_);

  template <class E>
  void operator+=(const PolynomialCRTP<E> &e_);

  template <class E>
  void operator-=(const PolynomialCRTP<E> &e_);

  void operator*=(double alpha);
  void operator/=(double alpha);

  template <class E>
  void deep_copy(const PolynomialCRTP<E> &e_);

  Cartesian<NVARS> operator()(const XY &xy) const;

  double a(int i, int j, int_t k) const;
  double &a(int i, int j, int_t k);

  double a(int_t i) const;
  double c(int_t i) const;

  const XY &x_center() const;
  double reference_length() const;

  double *coeffs_ptr();

protected:
  template <int DEG, int VARS>
  friend std::ostream &operator<<(std::ostream &os,
                                  const Poly2D<DEG, VARS> &poly2d);

private:
  double coeffs[n_vars() * n_coeffs()];
  double moments[n_coeffs()];

  int degree_;

  XY x_center_;
  double reference_length_;
};

/// Value that represents the smoothness of the polynomial.
template <int MAX_DEGREE, int NVARS>
Cartesian<NVARS> smoothness_indicator(const Poly2D<MAX_DEGREE, NVARS> &p);

template <int MAX_DEGREE, int NVARS>
std::ostream &operator<<(std::ostream &os,
                         const Poly2D<MAX_DEGREE, NVARS> &poly2d);

} // namespace zisa
#endif /* end of include guard */
