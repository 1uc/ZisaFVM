#ifndef GRID_DECL_H_8IQQ7
#define GRID_DECL_H_8IQQ7

#include <string>
#include <utility>

#include <zisa/math/cartesian.hpp>
#include <zisa/math/edge.hpp>
#include <zisa/math/triangle.hpp>
#include <zisa/memory/array.hpp>

namespace zisa {

struct Grid {
  int_t n_cells;
  int_t n_vertices;
  int_t n_edges;
  int_t n_interior_edges;
  int_t n_exterior_edges;

  int_t max_neighbours;

  array<int_t, 2> vertex_indices;
  array<int_t, 2> edge_indices;
  array<std::pair<int_t, int_t>, 1> left_right; ///< cell left & right of a edge

  array<int_t, 2> neighbours;
  array<bool, 2> is_valid;

  array<XY, 1> vertices;
  array<XY, 1> cell_centers;

  array<double, 1> volumes;
  array<XY, 1> normals;
  array<XY, 1> tangentials;

  array<array<double, 1>, 1> normalized_moments;

  Grid(array<XY, 1> vertices, array<int_t, 2> vertex_indices);

  const XY &vertex(int_t i, int_t k) const;
  Triangle triangle(int_t i) const;
  Edge edge(int_t e) const;
  Edge edge(int_t i, int_t k) const;
  double characteristic_length(int_t i) const;

  std::string str() const;
};

double volume(const Grid &grid);

double largest_circum_radius(const Grid &grid);

std::shared_ptr<Grid> load_gmsh(const std::string &filename);

/// Generate all moment for a 2D poly of degree 'deg'.
array<double, 1>
normalized_moments(const Triangle &tri, int deg, int_t quad_deg);

} // namespace zisa

#endif /* end of include guard */
