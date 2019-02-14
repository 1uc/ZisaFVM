#ifndef EDGE_RANGE_H_QYT12
#define EDGE_RANGE_H_QYT12

#include <zisa/grid/grid_decl.hpp>

namespace zisa {
class EdgeRange {
private:
  class EndIterator {};

  class Iterator {
  public:
    inline Iterator(const Grid &grid, int_t start, int_t end)
        : grid_(grid), current_(start), end_(end) {}

    inline void operator++() { ++current_; }

    inline std::pair<int_t, Edge> operator*() const {
      return {current_, grid_.edge(current_)};
    }

    inline bool operator!=(const EndIterator &) const {
      return current_ < end_;
    }
    inline bool operator==(const EndIterator &it) const {
      return !(*this != it);
    }

  private:
    const Grid &grid_;
    int_t current_;
    int_t end_;
  };

public:
  inline EdgeRange(const Grid &grid, int_t start, int_t end)
      : grid_(grid), start_(start), end_(end) {}

  inline Iterator begin() const { return Iterator(grid_, start_, end_); }
  inline EndIterator end() const { return EndIterator(); }

private:
  const Grid &grid_;
  int_t start_;
  int_t end_;
};

inline EdgeRange exterior_edges(const Grid &grid) {
  return EdgeRange(grid, grid.n_interior_edges, grid.n_edges);
}

inline EdgeRange interior_edges(const Grid &grid) {
  return EdgeRange(grid, 0, grid.n_interior_edges);
}

} // namespace zisa

#endif /* end of include guard */
