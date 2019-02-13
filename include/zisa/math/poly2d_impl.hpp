#ifndef POLY2D_IMPL_H_PP21T
#define POLY2D_IMPL_H_PP21T

#include "poly2d_decl.hpp"

namespace zisa {

constexpr ANY_DEVICE_INLINE int_t poly_dof(int deg) {
  return int_t((deg + 1) * (deg + 2) / 2);
}

constexpr ANY_DEVICE_INLINE int_t poly_index(int a, int b) {
  auto n = a + b;
  return int_t(((n + 1) * n) / 2 + b);
}

template <int MAX_DEGREE, int NVARS>
constexpr int Poly2D<MAX_DEGREE, NVARS>::max_degree() {
  return MAX_DEGREE;
}

template <int MAX_DEGREE, int NVARS>
int Poly2D<MAX_DEGREE, NVARS>::degree() const {
  return degree_;
}

template <int MAX_DEGREE, int NVARS>
constexpr int_t Poly2D<MAX_DEGREE, NVARS>::idx(int k, int l) {
  return poly_index(k, l);
}

template <int MAX_DEGREE, int NVARS>
constexpr int_t Poly2D<MAX_DEGREE, NVARS>::dof(int deg) {
  return poly_dof(deg);
}

template <int MAX_DEGREE, int NVARS>
constexpr int_t Poly2D<MAX_DEGREE, NVARS>::n_coeffs() {
  return dof(max_degree());
}

template <int MAX_DEGREE, int NVARS>
constexpr int_t Poly2D<MAX_DEGREE, NVARS>::n_vars() {
  return NVARS;
}

template <int MAX_DEGREE, int NVARS>
Poly2D<MAX_DEGREE, NVARS>::Poly2D() : degree_(0) {
  std::fill(coeffs, coeffs + n_coeffs() * n_vars(), 0.0);
  std::fill(moments, moments + n_coeffs(), 0.0);
  x_center_ = XY{0.0, 0.0};
  reference_length_ = 1.0;
}

template <int MAX_DEGREE, int NVARS>
Poly2D<MAX_DEGREE, NVARS>::Poly2D(int degree,
                                  const array<double, 1> &moments,
                                  const XY &x_center,
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

template <int MAX_DEGREE, int NVARS>
Poly2D<MAX_DEGREE, NVARS>::Poly2D(
    int degree,
    const std::initializer_list<double> &moments_list,
    const XY &x_center,
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

template <int MAX_DEGREE, int NVARS>
Poly2D<MAX_DEGREE, NVARS>::Poly2D(
    const std::initializer_list<double> &coeffs_list,
    const std::initializer_list<double> &moments_list,
    const XY &x_center,
    double reference_length)
    : degree_(poly_degree(coeffs_list.size())),
      x_center_(x_center),
      reference_length_(reference_length) {

  assert(n_vars() == 1);
  assert(coeffs_list.size() <= n_coeffs());
  assert(coeffs_list.size() == moments_list.size());
  assert(poly_dof(degree()) == coeffs_list.size());

  std::fill(coeffs, coeffs + n_vars() * n_coeffs(), 0.0);
  std::copy(coeffs_list.begin(), coeffs_list.end(), coeffs);

  std::fill(moments, moments + n_coeffs(), 0.0);
  std::copy(moments_list.begin(), moments_list.end(), moments);
}

template <int MAX_DEGREE, int NVARS>
template <class E>
Poly2D<MAX_DEGREE, NVARS>::Poly2D(const PolynomialCRTP<E> &e_) {
  deep_copy(static_cast<const E &>(e_));
}

template <int MAX_DEGREE, int NVARS>
template <class E>
void Poly2D<MAX_DEGREE, NVARS>::operator=(const PolynomialCRTP<E> &e_) {
  const E &e = static_cast<const E &>(e_);

  deep_copy(e);
}

template <int MAX_DEGREE, int NVARS>
template <class E>
void Poly2D<MAX_DEGREE, NVARS>::operator+=(const PolynomialCRTP<E> &e_) {
  const E &e = static_cast<const E &>(e_);

  *this = *this + e;
}

template <int MAX_DEGREE, int NVARS>
template <class E>
void Poly2D<MAX_DEGREE, NVARS>::operator-=(const PolynomialCRTP<E> &e_) {
  const E &e = static_cast<const E &>(e_);

  *this = *this - e;
}

template <int MAX_DEGREE, int NVARS>
void Poly2D<MAX_DEGREE, NVARS>::operator*=(double alpha) {
  *this = alpha * (*this);
}

template <int MAX_DEGREE, int NVARS>
void Poly2D<MAX_DEGREE, NVARS>::operator/=(double alpha) {
  *this *= 1.0 / alpha;
}

template <int MAX_DEGREE, int NVARS>
double *Poly2D<MAX_DEGREE, NVARS>::coeffs_ptr() {
  return coeffs;
}

template <int MAX_DEGREE, int NVARS>
template <class E>
void Poly2D<MAX_DEGREE, NVARS>::deep_copy(const PolynomialCRTP<E> &e_) {
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

template <int MAX_DEGREE, int NVARS>
Cartesian<NVARS> Poly2D<MAX_DEGREE, NVARS>::operator()(const XY &xy) const {
  auto x = (xy[0] - x_center_[0]) / reference_length_;
  auto y = (xy[1] - x_center_[1]) / reference_length_;

  auto d = degree();

  Cartesian<NVARS> px(0.0);
  double pow_x_kx = 1.0;
  for (int kx = 0; kx <= d; ++kx) {
    double pow_y_ky = 1.0;
    for (int ky = 0; ky <= d - kx; ++ky) {
      int_t i = idx(kx, ky);
      for (int_t k = 0; k < NVARS; ++k) {
        px[k] += a(i * NVARS + k) * (pow_x_kx * pow_y_ky - c(i));
      }
      pow_y_ky *= y;
    }

    pow_x_kx *= x;
  }

  return px;
}

template <int MAX_DEGREE, int NVARS>
double Poly2D<MAX_DEGREE, NVARS>::a(int i, int j, int_t k) const {
  return coeffs[idx(i, j) * NVARS + k];
}

template <int MAX_DEGREE, int NVARS>
double &Poly2D<MAX_DEGREE, NVARS>::a(int i, int j, int_t k) {
  return coeffs[idx(i, j) * NVARS + k];
}

template <int MAX_DEGREE, int NVARS>
double Poly2D<MAX_DEGREE, NVARS>::a(int_t i) const {
  return coeffs[i];
}

template <int MAX_DEGREE, int NVARS>
double Poly2D<MAX_DEGREE, NVARS>::c(int_t i) const {
  return moments[i];
}

template <int MAX_DEGREE, int NVARS>
const XY &Poly2D<MAX_DEGREE, NVARS>::x_center() const {
  return x_center_;
}

template <int MAX_DEGREE, int NVARS>
double Poly2D<MAX_DEGREE, NVARS>::reference_length() const {
  return reference_length_;
}

template <int MAX_DEGREE, int NVARS>
Cartesian<NVARS> smoothness_indicator(const Poly2D<MAX_DEGREE, NVARS> &p) {
  auto beta = Cartesian<NVARS>(0.0);

  int_t n = poly_dof(p.degree());
  for (int_t i = poly_dof(0); i < n; ++i) {
    for (int_t k = 0; k < p.n_vars(); ++k) {
      beta[k] += zisa::pow<2>(p.a(i * NVARS + k));
    }
  }

  return beta;
}

template <int MAX_DEGREE, int NVARS>
std::ostream &operator<<(std::ostream &os,
                         const Poly2D<MAX_DEGREE, NVARS> &poly2d) {
  os << format_as_list(poly2d.coeffs, poly_dof(poly2d.degree())) << "\n";
  os << format_as_list(poly2d.moments, poly_dof(poly2d.degree())) << "\n";

  return os;
}

} // namespace zisa
#endif /* end of include guard */
