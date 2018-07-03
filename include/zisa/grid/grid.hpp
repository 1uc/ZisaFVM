/* Unstructured 2D triangular grid.
 */

#ifndef GRID_H_RDTNK
#define GRID_H_RDTNK

#include <zisa/memory/array.hpp>
#include <zisa/math/cartesian.hpp>

namespace zisa {

struct Grid {
  array<int_t, 2> vertex_indices;
  array<int_t, 2> edge_indices;
  array<int_t, 2> neighbours;
  array<bool, 2> is_valid;

  array<XY, 1> vertices;
  array<XY, 1> cell_centers;

  array<double, 1> volumes;
  array<XY, 1> normals;
  array<XY, 1> tangentials;

  Grid(array<XY, 1> vertices, array<int_t, 2> vertex_indices);
};

std::shared_ptr<Grid> load_gmsh(const std::string &filename);

} // namespace zisa
#endif /* end of include guard */
