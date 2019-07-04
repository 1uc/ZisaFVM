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

template <class Derived, int MAX_DEGREE, int NVARS, int NDIMS>
class PolyND : public PolynomialCRTP<Derived> {
public:
  static constexpr int max_degree();
  int degree() const;

  static constexpr int_t dof(int deg);
  static constexpr int_t n_coeffs();
  static constexpr int_t n_vars();
  static constexpr int n_dims();

public:
  PolyND();

  PolyND(int degree,
         const array<double, 1> &moments,
         const XYZ &x_center,
         double reference_length);

  PolyND(int degree,
         const std::initializer_list<double> &moments_list,
         const XYZ &x_center,
         double reference_length);

  PolyND(const std::initializer_list<double> &coeffs_list,
         const std::initializer_list<double> &moments_list,
         const XYZ &x_center,
         double reference_length);

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

  template <class E>
  void deep_copy(const PolynomialCRTP<E> &e_);

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

  XYZ x_center_;
  double reference_length_;

private:
  template <class D, int MD, int NV, int ND>
  friend std::ostream &operator<<(std::ostream &os,
                                  const PolyND<D, MD, NV, ND> &poly);
};

template <class Derived, int MAX_DEGREE, int NVARS, int NDIMS>
std::ostream &operator<<(std::ostream &os,
                         const PolyND<Derived, MAX_DEGREE, NVARS, NDIMS> &poly);

template <int MAX_DEGREE, int NVARS>
class Poly2D : public PolyND<Poly2D<MAX_DEGREE, NVARS>, MAX_DEGREE, NVARS, 2> {
private:
  using super = PolyND<Poly2D<MAX_DEGREE, NVARS>, MAX_DEGREE, NVARS, 2>;

public:
  static constexpr int_t idx(int k, int l);

public:
  using super::super;
  using super::operator=;
  using super::operator+=;
  using super::operator-=;
  using super::operator*=;
  using super::operator/=;

  Cartesian<NVARS> operator()(const XYZ &xy) const;
};

template <int MAX_DEGREE, int NVARS>
class Poly3D : public PolyND<Poly3D<MAX_DEGREE, NVARS>, MAX_DEGREE, NVARS, 3> {
private:
  using super = PolyND<Poly3D<MAX_DEGREE, NVARS>, MAX_DEGREE, NVARS, 3>;

public:
  static constexpr int_t idx(int i, int j, int k);

public:
  using super::super;
  using super::operator=;
  using super::operator+=;
  using super::operator-=;
  using super::operator*=;
  using super::operator/=;

  Cartesian<NVARS> operator()(const XYZ &xy) const;
};

/// Value that represents the smoothness of the polynomial.
template <class Derived, int MAX_DEGREE, int NVARS, int NDIMS>
Cartesian<NVARS>
smoothness_indicator(const PolyND<Derived, MAX_DEGREE, NVARS, NDIMS> &p);

} // namespace zisa
#endif /* end of include guard */
