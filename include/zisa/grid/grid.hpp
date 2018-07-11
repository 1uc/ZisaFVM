/* Unstructured 2D triangular grid.
 */

#ifndef GRID_H_RDTNK
#define GRID_H_RDTNK

#include <zisa/math/cartesian.hpp>
#include <zisa/math/triangle.hpp>
#include <zisa/memory/array.hpp>

namespace zisa {

struct Grid {
  int_t n_cells;
  int_t max_neighbours;

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

  Triangle triangles(int_t i) const;
};

class TriangleRange {
private:
  class EndIterator {};

  class Iterator {
  public:
    inline Iterator(const Grid &grid) : grid(grid), i(0) {}

    inline void operator++() { i++; }
    inline Triangle operator*() const { return grid.triangles(i); }

    inline bool operator!=(const EndIterator &) const {
      return i < grid.n_cells;
    }
    inline bool operator==(const EndIterator &it) const {
      return !(*this != it);
    }

  private:
    const Grid &grid;
    int_t i;
  };

public:
  inline TriangleRange(const Grid &grid) : grid(grid) {}

  inline Iterator begin() const { return Iterator(grid); }
  inline EndIterator end() const { return EndIterator(); }

private:
  const Grid &grid;
};

inline TriangleRange triangles(const Grid &grid) { return TriangleRange(grid); }

double largest_circum_radius(const Grid &grid);

std::shared_ptr<Grid> load_gmsh(const std::string &filename);

} // namespace zisa
#endif /* end of include guard */
