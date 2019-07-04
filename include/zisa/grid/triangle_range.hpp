#ifndef TRIANGLE_RANGE_H_51FI7
#define TRIANGLE_RANGE_H_51FI7

#include <utility>

#include <zisa/grid/grid_decl.hpp>

namespace zisa {

class TriangleRange {
private:
  class EndIterator {};

  class Iterator {
  public:
    explicit inline Iterator(const Grid &grid) : grid(grid), i(0) {}

    inline void operator++() { i++; }
    inline std::pair<int_t, Triangle> operator*() const {
      return {i, triangle(grid, i)};
    }

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
  explicit inline TriangleRange(const Grid &grid) : grid(grid) {}

  static constexpr bool has_item() { return true; }
  inline Triangle item(int_t i) const { return triangle(grid, i); }

  inline Iterator begin() const { return Iterator(grid); }
  inline EndIterator end() const { return EndIterator(); }

  inline int_t start_index() const { return 0; }
  inline int_t end_index() const { return grid.n_cells; }

private:
  const Grid &grid;
};

inline TriangleRange triangles(const Grid &grid) { return TriangleRange(grid); }

} // namespace zisa
#endif /* end of include guard */
