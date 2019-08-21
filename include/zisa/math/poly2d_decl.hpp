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
template <int NDIMS>
constexpr ANY_DEVICE_INLINE int_t poly_dof(int deg);

/// Highest degree of a polynomial with '<=n_coeffs' coefficients.
template <int NDIMS>
int poly_degree(int_t n_coeffs);

/// Highest degree of a polynomial with '<=n_coeffs' coefficients.
int poly_degree(int_t n_coeffs, int n_dims);

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
 *    @param a  -- slow index
 *    @param b  -- fast index, i.e. for fixed (k+l) consecutive values of l
 *                 result in consecutive values of the linear index.
 **/
constexpr ANY_DEVICE_INLINE int_t poly_index(int a, int b);

/// Compute linear index of the coefficient of a 3D polynomial.
constexpr ANY_DEVICE_INLINE int_t poly_index(int a, int b, int c);

template <int NCOEFFS, int NVARS>
class PolyND : public PolynomialCRTP<PolyND<NCOEFFS, NVARS>> {
public:
  int max_degree() const;
  int degree() const;

  static constexpr int_t n_coeffs();
  static constexpr int_t n_vars();

  static constexpr int_t idx(int k, int l);
  static constexpr int_t idx(int k, int l, int m);

public:
  PolyND();

  PolyND(int degree,
         const array<double, 1> &moments,
         const XYZ &x_center,
         double reference_length,
         int n_dims);

  PolyND(int degree,
         const std::initializer_list<double> &moments_list,
         const XYZ &x_center,
         double reference_length,
         int n_dims);

  PolyND(const std::initializer_list<double> &coeffs_list,
         const std::initializer_list<double> &moments_list,
         const XYZ &x_center,
         double reference_length,
         int n_dims);

  template <class E>
  PolyND(const PolynomialCRTP<E> &e_);

  template <class E>
  void operator=(const PolynomialCRTP<E> &e_);

  template <class E>
  void operator+=(const PolynomialCRTP<E> &e_);

  template <class E>
  void operator-=(const PolynomialCRTP<E> &e_);

  void operator*=(double alpha);
  void operator/=(double alpha);

  Cartesian<NVARS> operator()(const XYZ &xyz) const;

  Cartesian<NVARS> eval_2d(const XYZ &xyz) const;
  Cartesian<NVARS> eval_3d(const XYZ &xyz) const;

  template <class E>
  void deep_copy(const PolynomialCRTP<E> &e_);

  int n_dims() const;
  int_t dof() const;

  double a(int_t i) const;
  double &a(int_t i);

  double c(int_t i) const;

  const XYZ &x_center() const;
  double reference_length() const;

  double *coeffs_ptr();

private:
  double coeffs[n_vars() * n_coeffs()];
  double moments[n_coeffs()];

  int degree_;
  int n_dims_;

  XYZ x_center_;
  double reference_length_;

private:
  template <int M, int NV>
  friend std::ostream &operator<<(std::ostream &os, const PolyND<M, NV> &poly);
};

template <int NCOEFFS, int NVARS>
std::ostream &operator<<(std::ostream &os, const PolyND<NCOEFFS, NVARS> &poly);

/// Value that represents the smoothness of the polynomial.
template <int NCOEFFS, int NVARS>
Cartesian<NVARS> smoothness_indicator(const PolyND<NCOEFFS, NVARS> &p);

} // namespace zisa
#endif /* end of include guard */
