#ifndef GRID_DECL_H_8IQQ7
#define GRID_DECL_H_8IQQ7

#include <optional>
#include <string>
#include <utility>

#include <zisa/grid/gmsh_reader.hpp>
#include <zisa/io/hdf5_writer_fwd.hpp>
#include <zisa/math/cartesian.hpp>
#include <zisa/math/cell.hpp>
#include <zisa/math/denormalized_rule.hpp>
#include <zisa/math/edge.hpp>
#include <zisa/math/face.hpp>
#include <zisa/math/tetrahedron.hpp>
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

  array<XYZ, 1> vertices;
  array<XYZ, 1> cell_centers;
  array<XYZ, 1> face_centers;

  array<Cell, 1> cells;
  array<Face, 1> faces;

  array<double, 1> volumes;
  array<double, 1> inradii;
  array<double, 1> circum_radii;
  array<XYZ, 1> normals;
  array<XYZ, 2> tangentials;

  array<array<double, 1>, 1> normalized_moments;

  Grid() = default;

  /// Generate a grid optionally with quadrature rules.
  Grid(GMSHElementType element_type,
       array<XYZ, 1> vertices,
       array<int_t, 2> vertex_indices,
       int_t quad_deg = 0);

  [[nodiscard]] static Grid load(HDF5Reader &reader);

  const XYZ &vertex(int_t i, int_t k) const;
  Triangle triangle(int_t i) const;
  Edge edge(int_t e) const;
  Edge edge(int_t i, int_t k) const;

  Face face(int_t i, int_t k) const;
  const XYZ &face_center(int_t i, int_t k) const;

  double inradius(int_t i) const;
  double circum_radius(int_t i) const;
  double characteristic_length(int_t i) const;

  std::string str() const;
};

void save(HDF5Writer &writer, const Grid &grid);

double volume(const Grid &grid);
Triangle triangle(const Grid &grid, int_t i);

bool is_inside_cell(const Grid &grid, int_t i, const XYZ &x);

std::optional<int_t> locate(const Grid &grid, const XYZ &x, int_t i_guess = 0);

double largest_circum_radius(const Grid &grid);
double smallest_inradius(const Grid &grid);

std::shared_ptr<Grid> load_grid(const std::string &filename, int_t quad_deg);
std::shared_ptr<Grid> load_gmsh(const std::string &filename, int_t quad_deg);
std::shared_ptr<Grid> load_gmsh(const std::string &filename);

/// Generate all moment for a 2D poly of degree 'deg'.
array<double, 1>
normalized_moments(const Triangle &tri, int deg, int_t quad_deg);

/// Generate all moment for a 3D poly of degree 'deg'.
array<double, 1>
normalized_moments(const Tetrahedron &tet, int deg, int_t quad_deg);

} // namespace zisa

#endif /* end of include guard */
