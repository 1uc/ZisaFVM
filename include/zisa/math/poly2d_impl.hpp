#ifndef POLY2D_IMPL_H_PP21T
#define POLY2D_IMPL_H_PP21T

#include "poly2d_decl.hpp"

namespace zisa {

template <int NDIMS>
constexpr ANY_DEVICE_INLINE int_t poly_dof(int deg) {
  if constexpr (NDIMS == 2) {
    return int_t(((deg + 1) * (deg + 2)) / 2);
  } else {
    return int_t(((deg + 1) * (deg + 2) * (deg + 3)) / 6);
  }
}

template <int NDIMS>
int poly_degree(int_t n_coeffs) {

  LOG_ERR_IF(n_coeffs == 0, "Too few coefficients.");

  for (int d = 1;; ++d) {
    if (poly_dof<NDIMS>(d) > n_coeffs) {
      return d - 1;
    }
  }
}

constexpr ANY_DEVICE_INLINE int_t poly_index(int a, int b) {
  auto n = a + b;
  return int_t(((n + 1) * n) / 2 + b);
}

constexpr ANY_DEVICE_INLINE int_t poly_index(int a, int b, int c) {
  return poly_dof<3>(a + b + c - 1) + poly_index(b, c);
}

template <class Derived, int MAX_DEGREE, int NVARS, int NDIMS>
constexpr int PolyND<Derived, MAX_DEGREE, NVARS, NDIMS>::max_degree() {
  return MAX_DEGREE;
}

template <class Derived, int MAX_DEGREE, int NVARS, int NDIMS>
int PolyND<Derived, MAX_DEGREE, NVARS, NDIMS>::degree() const {
  return degree_;
}

template <class Derived, int MAX_DEGREE, int NVARS, int NDIMS>
constexpr int_t PolyND<Derived, MAX_DEGREE, NVARS, NDIMS>::dof(int deg) {
  return poly_dof<NDIMS>(deg);
}

template <class Derived, int MAX_DEGREE, int NVARS, int NDIMS>
constexpr int_t PolyND<Derived, MAX_DEGREE, NVARS, NDIMS>::n_coeffs() {
  return dof(max_degree());
}

template <class Derived, int MAX_DEGREE, int NVARS, int NDIMS>
constexpr int_t PolyND<Derived, MAX_DEGREE, NVARS, NDIMS>::n_vars() {
  return NVARS;
}

template <class Derived, int MAX_DEGREE, int NVARS, int NDIMS>
constexpr int PolyND<Derived, MAX_DEGREE, NVARS, NDIMS>::n_dims() {
  return NDIMS;
}

template <class Derived, int MAX_DEGREE, int NVARS, int NDIMS>
PolyND<Derived, MAX_DEGREE, NVARS, NDIMS>::PolyND() : degree_(0) {
  std::fill(coeffs, coeffs + n_coeffs() * n_vars(), 0.0);
  std::fill(moments, moments + n_coeffs(), 0.0);
  x_center_ = XYZ::zeros();
  reference_length_ = 1.0;
}

template <class Derived, int MAX_DEGREE, int NVARS, int NDIMS>
PolyND<Derived, MAX_DEGREE, NVARS, NDIMS>::PolyND(
    int degree,
    const array<double, 1> &moments,
    const XYZ &x_center,
    double reference_length)
    : degree_(degree),
      x_center_(x_center),
      reference_length_(reference_length) {

  std::fill(this->coeffs, this->coeffs + n_coeffs() * n_vars(), 0.0);

  static_assert(n_coeffs() >= 3, "Unusually low degree; fix code below.");
  this->moments[0] = this->moments[1] = this->moments[2] = 0.0;

  if (std::distance(moments.begin(), moments.end()) > 3) {
    std::copy(moments.begin() + 3, moments.end(), this->moments + 3);
  }
}

template <class Derived, int MAX_DEGREE, int NVARS, int NDIMS>
PolyND<Derived, MAX_DEGREE, NVARS, NDIMS>::PolyND(
    int degree,
    const std::initializer_list<double> &moments_list,
    const XYZ &x_center,
    double reference_length)
    : degree_(degree),
      x_center_(x_center),
      reference_length_(reference_length) {

  std::fill(this->coeffs, this->coeffs + n_coeffs() * n_vars(), 0.0);

  static_assert(n_coeffs() >= 3, "Unusually low degree; fix code below.");
  this->moments[0] = this->moments[1] = this->moments[2] = 0.0;

  if (std::distance(moments_list.begin(), moments_list.end()) > 3) {
    std::copy(moments_list.begin() + 3, moments_list.end(), moments + 3);
  }
}

template <class Derived, int MAX_DEGREE, int NVARS, int NDIMS>
PolyND<Derived, MAX_DEGREE, NVARS, NDIMS>::PolyND(
    const std::initializer_list<double> &coeffs_list,
    const std::initializer_list<double> &moments_list,
    const XYZ &x_center,
    double reference_length)
    : degree_(poly_degree<NDIMS>(coeffs_list.size())),
      x_center_(x_center),
      reference_length_(reference_length) {

  assert(n_vars() == 1);
  assert(coeffs_list.size() <= n_coeffs());
  assert(coeffs_list.size() == moments_list.size());
  assert(poly_dof<NDIMS>(degree()) == coeffs_list.size());

  std::fill(coeffs, coeffs + n_vars() * n_coeffs(), 0.0);
  std::copy(coeffs_list.begin(), coeffs_list.end(), coeffs);

  std::fill(moments, moments + n_coeffs(), 0.0);
  std::copy(moments_list.begin(), moments_list.end(), moments);
}

template <class Derived, int MAX_DEGREE, int NVARS, int NDIMS>
template <class E>
PolyND<Derived, MAX_DEGREE, NVARS, NDIMS>::PolyND(const PolynomialCRTP<E> &e_) {
  const auto &e = static_cast<const E &>(e_);

  static_assert(E::n_dims() == NDIMS, "Dimensions mismatch.");
  deep_copy(e);
}

template <class Derived, int MAX_DEGREE, int NVARS, int NDIMS>
double *PolyND<Derived, MAX_DEGREE, NVARS, NDIMS>::coeffs_ptr() {
  return coeffs;
}

template <class Derived, int MAX_DEGREE, int NVARS, int NDIMS>
template <class E>
void PolyND<Derived, MAX_DEGREE, NVARS, NDIMS>::deep_copy(
    const PolynomialCRTP<E> &e_) {
  const E &e = static_cast<const E &>(e_);

  static_assert(E::max_degree() <= max_degree(), "Polynomial degree mismatch.");

  degree_ = e.degree();

  constexpr int_t n = n_coeffs() * n_vars();
  for (int_t i = 0; i < n; ++i) {
    coeffs[i] = e.a(i);
  }

  for (int_t i = 0; i < n_coeffs(); ++i) {
    moments[i] = e.c(i);
  }

  x_center_ = e.x_center();
  reference_length_ = e.reference_length();
}

template <class Derived, int MAX_DEGREE, int NVARS, int NDIMS>
double PolyND<Derived, MAX_DEGREE, NVARS, NDIMS>::a(int_t i) const {
  return coeffs[i];
}

template <class Derived, int MAX_DEGREE, int NVARS, int NDIMS>
double &PolyND<Derived, MAX_DEGREE, NVARS, NDIMS>::a(int_t i) {
  return coeffs[i];
}

template <class Derived, int MAX_DEGREE, int NVARS, int NDIMS>
double PolyND<Derived, MAX_DEGREE, NVARS, NDIMS>::c(int_t i) const {
  return moments[i];
}

template <class Derived, int MAX_DEGREE, int NVARS, int NDIMS>
const XYZ &PolyND<Derived, MAX_DEGREE, NVARS, NDIMS>::x_center() const {
  return x_center_;
}

template <class Derived, int MAX_DEGREE, int NVARS, int NDIMS>
double PolyND<Derived, MAX_DEGREE, NVARS, NDIMS>::reference_length() const {
  return reference_length_;
}

template <class Derived, int MAX_DEGREE, int NVARS, int NDIMS>
template <class E>
void PolyND<Derived, MAX_DEGREE, NVARS, NDIMS>::
operator=(const PolynomialCRTP<E> &e_) {
  static_assert(E::n_dims() == NDIMS, "Dimensions mismatch.");

  const E &e = static_cast<const E &>(e_);
  this->deep_copy(e);
}

template <class Derived, int MAX_DEGREE, int NVARS, int NDIMS>
template <class E>
void PolyND<Derived, MAX_DEGREE, NVARS, NDIMS>::
operator+=(const PolynomialCRTP<E> &e_) {
  static_assert(E::n_dims() == NDIMS, "Dimensions mismatch.");

  const E &e = static_cast<const E &>(e_);
  *this = *this + e;
}

template <class Derived, int MAX_DEGREE, int NVARS, int NDIMS>
template <class E>
void PolyND<Derived, MAX_DEGREE, NVARS, NDIMS>::
operator-=(const PolynomialCRTP<E> &e_) {
  static_assert(E::n_dims() == NDIMS, "Dimensions mismatch.");

  const E &e = static_cast<const E &>(e_);
  *this = *this - e;
}

template <class Derived, int MAX_DEGREE, int NVARS, int NDIMS>
void PolyND<Derived, MAX_DEGREE, NVARS, NDIMS>::operator*=(double alpha) {
  *this = alpha * (*this);
}

template <class Derived, int MAX_DEGREE, int NVARS, int NDIMS>
void PolyND<Derived, MAX_DEGREE, NVARS, NDIMS>::operator/=(double alpha) {
  *this *= 1.0 / alpha;
}

template <int MAX_DEGREE, int NVARS>
constexpr int_t Poly2D<MAX_DEGREE, NVARS>::idx(int k, int l) {
  return poly_index(k, l);
}

template <int MAX_DEGREE, int NVARS>
Cartesian<NVARS> Poly2D<MAX_DEGREE, NVARS>::operator()(const XYZ &xy) const {
  auto [x, y, z] = XYZ((xy - this->x_center()) / this->reference_length());

  auto d = this->degree();

  Cartesian<NVARS> px(0.0);
  double pow_x_kx = 1.0;
  for (int kx = 0; kx <= d; ++kx) {
    double pow_y_ky = 1.0;
    for (int ky = 0; ky <= d - kx; ++ky) {
      int_t i = idx(kx, ky);
      for (int_t k = 0; k < NVARS; ++k) {
        px[k] += this->a(i * NVARS + k) * (pow_x_kx * pow_y_ky - this->c(i));
      }
      pow_y_ky *= y;
    }

    pow_x_kx *= x;
  }

  return px;
}

template <int MAX_DEGREE, int NVARS>
constexpr int_t Poly3D<MAX_DEGREE, NVARS>::idx(int i, int j, int k) {
  return poly_index(i, j, k);
}

template <int MAX_DEGREE, int NVARS>
Cartesian<NVARS> Poly3D<MAX_DEGREE, NVARS>::operator()(const XYZ &xy) const {
  auto [x, y, z] = XYZ((xy - this->x_center()) / this->reference_length());

  auto d = this->degree();

  Cartesian<NVARS> px(0.0);
  double pow_x_kx = 1.0;
  for (int kx = 0; kx <= d; ++kx) {
    double pow_y_ky = 1.0;
    for (int ky = 0; ky <= d - kx; ++ky) {
      double pow_z_kz = 1.0;

      for (int kz = 0; kz <= d - kx - ky; ++kz) {
        int_t i = idx(kx, ky, kz);
        for (int_t k = 0; k < NVARS; ++k) {
          double ak = this->a(i * NVARS + k);
          px[k] += ak * (pow_x_kx * pow_y_ky * pow_z_kz - this->c(i));
        }

        pow_z_kz *= z;
      }
      pow_y_ky *= y;
    }
    pow_x_kx *= x;
  }

  return px;
}

template <class Derived, int MAX_DEGREE, int NVARS, int NDIMS>
Cartesian<NVARS>
smoothness_indicator(const PolyND<Derived, MAX_DEGREE, NVARS, NDIMS> &p) {
  auto beta = Cartesian<NVARS>(0.0);

  int_t n = poly_dof<NDIMS>(p.degree());
  for (int_t i = poly_dof<NDIMS>(0); i < n; ++i) {
    for (int_t k = 0; k < p.n_vars(); ++k) {
      beta[k] += zisa::pow<2>(p.a(i * NVARS + k));
    }
  }

  return beta;
}

template <class Derived, int MAX_DEGREE, int NVARS, int NDIMS>
std::ostream &
operator<<(std::ostream &os,
           const PolyND<Derived, MAX_DEGREE, NVARS, NDIMS> &poly) {
  os << format_as_list(poly.coeffs, poly_dof<NDIMS>(poly.degree())) << "\n";
  os << format_as_list(poly.moments, poly_dof<NDIMS>(poly.degree())) << "\n";

  return os;
}

} // namespace zisa
#endif /* end of include guard */
