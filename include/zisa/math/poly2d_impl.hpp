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

template <int MAX_DEGREE>
constexpr int Poly2D<MAX_DEGREE>::max_degree() {
  return MAX_DEGREE;
}

template <int MAX_DEGREE>
int Poly2D<MAX_DEGREE>::degree() const {
  return degree_;
}

template <int MAX_DEGREE>
constexpr int_t Poly2D<MAX_DEGREE>::idx(int k, int l) {
  return poly_index(k, l);
}

template <int MAX_DEGREE>
constexpr int_t Poly2D<MAX_DEGREE>::dof(int deg) {
  return poly_dof(deg);
}

template <int MAX_DEGREE>
constexpr int_t Poly2D<MAX_DEGREE>::n_coeffs() {
  return dof(max_degree());
}

template <int MAX_DEGREE>
Poly2D<MAX_DEGREE>::Poly2D() : degree_(0) {
  std::fill(coeffs, coeffs + n_coeffs(), 0.0);
  std::fill(moments, moments + n_coeffs(), 0.0);
  x_center_ = XY{0.0, 0.0};
  reference_length_ = 1.0;
}

template <int MAX_DEGREE>
Poly2D<MAX_DEGREE>::Poly2D(const std::initializer_list<double> &coeffs_list,
                           const std::initializer_list<double> &moments_list,
                           const XY &x_center,
                           double reference_length)
    : degree_(poly_degree(coeffs_list.size())),
      x_center_(x_center),
      reference_length_(reference_length) {

  assert(coeffs_list.size() <= n_coeffs());
  assert(coeffs_list.size() == moments_list.size());
  assert(poly_dof(degree()) == coeffs_list.size());

  std::fill(coeffs, coeffs + n_coeffs(), 0.0);
  std::copy(coeffs_list.begin(), coeffs_list.end(), coeffs);
  std::copy(moments_list.begin(), moments_list.end(), moments);
}

template <int MAX_DEGREE>
template <class E>
Poly2D<MAX_DEGREE>::Poly2D(const PolynomialCRTP<E> &e_) {
  deep_copy(static_cast<const E &>(e_));
}

template <int MAX_DEGREE>
template <int D>
void Poly2D<MAX_DEGREE>::operator=(const Poly2D<D> &other) {
  deep_copy(other);
}

template <int MAX_DEGREE>
template <class E>
void Poly2D<MAX_DEGREE>::operator=(const PolynomialCRTP<E> &e_) {
  const E &e = static_cast<const E &>(e_);

  deep_copy(e);
}

template <int MAX_DEGREE>
template <class E>
void Poly2D<MAX_DEGREE>::operator+=(const PolynomialCRTP<E> &e_) {
  const E &e = static_cast<const E &>(e_);

  *this = *this + e;
}

template <int MAX_DEGREE>
template <class E>
void Poly2D<MAX_DEGREE>::operator-=(const PolynomialCRTP<E> &e_) {
  const E &e = static_cast<const E &>(e_);

  *this = *this - e;
}

template <int MAX_DEGREE>
void Poly2D<MAX_DEGREE>::operator*=(double alpha) {
  *this = alpha * (*this);
}

template <int MAX_DEGREE>
void Poly2D<MAX_DEGREE>::operator/=(double alpha) {
  *this *= 1.0 / alpha;
}

template <int MAX_DEGREE>
template <class E>
void Poly2D<MAX_DEGREE>::deep_copy(const PolynomialCRTP<E> &e_) {
  const E &e = static_cast<const E &>(e_);

  static_assert(E::max_degree() <= max_degree(), "Polynomial degree mismatch.");

  degree_ = e.degree();

  int deg = degree();
  for (int d = 0; d <= deg; ++d) {
    for (int k = 0; k <= d; ++k) {
      int l = d - k;
      this->a(k, l) = e.a(k, l);
      this->c(k, l) = e.c(k, l);
    }
  }

  offset_cached = false;

  x_center_ = e.x_center();
  reference_length_ = e.reference_length();
}

template <class Coeffs>
double horners_method(const Coeffs &coeffs, double x, int n) {
  double px = coeffs(n);
  for (int k = n - 1; k >= 0; --k) {
    px = x * px + coeffs(k);
  }

  return px;
}

template <int MAX_DEGREE>
void Poly2D<MAX_DEGREE>::cache_offset() const {

  auto d = degree();

  double tmp = 0.0;
  for (int k = 0; k <= d; ++k) {
    for (int l = 0; l <= d - k; ++l) {
      tmp += a(k, l) * c(k, l);
    }
  }

  offset = tmp;
  offset_cached = true;
}

template <int MAX_DEGREE>
double Poly2D<MAX_DEGREE>::operator()(const XY &xy) const {
  auto x = (xy[0] - x_center_[0]) / reference_length_;
  auto y = (xy[1] - x_center_[1]) / reference_length_;

  auto d = degree();

  if (!offset_cached) {
    cache_offset();
  }

  auto px = horners_method(
      [this, d, y](int k) {
        return horners_method([this, k](int l) { return a(k, l); }, y, d - k);
      },
      x,
      d);

  return px - offset;
}

template <int MAX_DEGREE>
double &Poly2D<MAX_DEGREE>::a(int i, int j) {
  return coeffs[idx(i, j)];
}

template <int MAX_DEGREE>
double &Poly2D<MAX_DEGREE>::c(int i, int j) {
  return moments[idx(i, j)];
}

template <int MAX_DEGREE>
double Poly2D<MAX_DEGREE>::a(int i, int j) const {
  if (i + j <= degree()) {
    return coeffs[idx(i, j)];
  } else {
    return 0.0;
  }
}

template <int MAX_DEGREE>
double Poly2D<MAX_DEGREE>::c(int i, int j) const {
  return moments[idx(i, j)];
}

template <int MAX_DEGREE>
const XY &Poly2D<MAX_DEGREE>::x_center() const {
  return x_center_;
}

template <int MAX_DEGREE>
double Poly2D<MAX_DEGREE>::reference_length() const {
  return reference_length_;
}

template <int MAX_DEGREE>
double smoothness_indicator(const Poly2D<MAX_DEGREE> &p) {
  double beta = 0.0;

  for (int d = 1; d <= p.degree(); ++d) {
    for (int k = 0; k <= d; ++k) {
      beta += zisa::pow<2>(p.a(k, d - k));
    }
  }

  return beta;
}

template <int MAX_DEGREE>
std::ostream &operator<<(std::ostream &os, const Poly2D<MAX_DEGREE> &poly2d) {
  os << format_as_list(poly2d.coeffs, poly_dof(poly2d.degree())) << "\n";
  os << format_as_list(poly2d.moments, poly_dof(poly2d.degree())) << "\n";

  return os;
}

} // namespace zisa
#endif /* end of include guard */
