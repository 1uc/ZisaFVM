/*
 *
 */

#ifndef WENO_AO_H_A1UM2
#define WENO_AO_H_A1UM2

#include <zisa/config.hpp>
#include <zisa/grid/grid.hpp>
#include <zisa/math/cartesian.hpp>

namespace zisa {

struct LocalIndex {
  int_t local;
  int_t global;
};

class WENO_AO {
public:
  WENO_AO(const Grid &grid, int_t i_cell);

private:
  array<array<LocalIndex, 1>, 1> stencils;
  int max_order;
};

class Cone {
public:
  Cone(const XY &x, const XY &dir, double cos_angle);
  bool is_inside(const XY &y) const;

private:
  XY x;
  XY dir;
  double cos_angle;
};

std::vector<int_t> central_stencil(const Grid &grid, int_t i, int_t n_points);
std::vector<int_t>
biased_stencil(const Grid &grid, int_t i, int_t n_points, const Cone &cone);

} // namespace zisa
#endif /* end of include guard */
