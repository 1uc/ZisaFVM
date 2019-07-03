#ifndef ZISA_CELL_RANGE_HPP
#define ZISA_CELL_RANGE_HPP

#include <zisa/config.hpp>

#include <zisa/grid/grid.hpp>
#include <zisa/loops/range.hpp>
#include <zisa/math/cell.hpp>

namespace zisa {

class CellRange {
private:
  class EndIterator {};

  class Iterator {
  public:
    explicit inline Iterator(const Grid &grid) : grid(grid), i(0) {}

    inline void operator++() { i++; }
    inline std::pair<int_t, Cell> operator*() const {
      return std::pair<int_t, Cell>{i, grid.cells(i)};
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
  explicit inline CellRange(const Grid &grid) : grid(grid) {}

  static constexpr bool has_item() { return true; }
  inline Cell item(int_t i) const { return grid.cells(i); }

  inline Iterator begin() const { return Iterator(grid); }
  inline EndIterator end() const { return EndIterator(); }

  inline int_t start_index() const { return 0; }
  inline int_t end_index() const { return grid.n_cells; }

private:
  const Grid &grid;
};

inline CellRange cells(const Grid &grid) { return CellRange(grid); }

inline PlainIndexRange cell_indices(const Grid &grid) {
  return PlainIndexRange(0, grid.n_cells);
}

inline PlainIndexRange interior_face_indices(const Grid &grid) {
  return PlainIndexRange(0, grid.n_interior_edges);
}

}
#endif // ZISA_CELL_RANGE_HPP
