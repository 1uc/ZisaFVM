/*
 *
 */

#include <zisa/math/poly2d.hpp>

namespace zisa {

ANY_DEVICE int_t poly_dof(int deg, int n_dims) {
  if (n_dims == 2) {
    return poly_dof<2>(deg);
  } else {
    return poly_dof<3>(deg);
  }
}

} // namespace zisa
