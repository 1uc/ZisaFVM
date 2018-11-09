#ifndef GRID_VARIABLES_IMPL_H_FYDQZ
#define GRID_VARIABLES_IMPL_H_FYDQZ

#include "grid_variables_decl.hpp"

#include <zisa/math/cartesian_expr.hpp>

namespace zisa {

class CVarExpr : public CartesianExpr<CVarExpr, double> {
public:
  CVarExpr(GridVariables &array, int_t i) : array(array), i(i) {}

  double &operator[](int_t k) { return array(i, k); }
  double operator()(int_t k) const { return array(i, k); }

  template <class E>
  void operator=(const CartesianExpr<E, double> &e_) {
    const E &e = static_cast<const E &>(e_);

    for (int_t k = 0; k < E::size(); ++k) {
      array(i, k) = e(k);
    }
  }

  int_t size() const { return array.shape(1); }

private:
  GridVariables &array;
  int_t i;
};

class CVarConstExpr : public CartesianExpr<CVarConstExpr, double> {
public:
  CVarConstExpr(const GridVariables &array, int_t i) : array(array), i(i) {}
  double operator[](int_t k) const { return array(i, k); }
  double operator()(int_t k) const { return array(i, k); }

  int_t size() const { return array.shape(1); }

private:
  const GridVariables &array;
  int_t i;
};

inline CVarExpr GridVariables::operator()(int_t i) {
  return CVarExpr(*this, i);
}
inline CVarConstExpr GridVariables::operator()(int_t i) const {
  return CVarConstExpr(*this, i);
}

inline double &GridVariables::operator()(int_t i, int_t k) {
  return super::operator()(i, k);
}

inline double GridVariables::operator()(int_t i, int_t k) const {
  return super::operator()(i, k);
}

} // namespace zisa
#endif /* end of include guard */
