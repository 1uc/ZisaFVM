#ifndef EDGE_RANGE_H_QYT12
#define EDGE_RANGE_H_QYT12

#include <zisa/grid/grid_decl.hpp>

namespace zisa {

class EdgeRange {
private:
  class EndIterator {};

  class Iterator {
  public:
    inline Iterator(const Grid &grid) : grid(grid), e(0) {}

    inline void operator++() { ++e; }

    inline std::pair<int_t, Edge> operator*() const {
      return {e, grid.edge(e)};
    }

    inline bool operator!=(const EndIterator &) const {
      return e < grid.n_edges;
    }
    inline bool operator==(const EndIterator &it) const {
      return !(*this != it);
    }

  private:
    const Grid &grid;
    int_t e;
  };

public:
  inline EdgeRange(const Grid &grid) : grid(grid) {}

  inline Iterator begin() const { return Iterator(grid); }
  inline EndIterator end() const { return EndIterator(); }

private:
  const Grid &grid;
};

inline EdgeRange interior_edges(const Grid &grid) { return EdgeRange(grid); }

} // namespace zisa

#endif /* end of include guard */
