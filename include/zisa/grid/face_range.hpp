// SPDX-License-Identifier: MIT
// Copyright (c) 2021 ETH Zurich, Luc Grosheintz-Laval

#ifndef ZISA_FACE_RANGE_HPP_NUCPE
#define ZISA_FACE_RANGE_HPP_NUCPE

#include <zisa/config.hpp>
#include <zisa/grid/grid_decl.hpp>

namespace zisa {
class FaceRange {
private:
  class EndIterator {};

  class Iterator {
  public:
    inline Iterator(const Grid &grid, int_t start, int_t end)
        : grid_(grid), current_(start), end_(end) {}

    inline void operator++() { ++current_; }

    inline std::pair<int_t, Face> operator*() const {
      return {current_, grid_.faces(current_)};
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
  inline FaceRange(const Grid &grid, int_t start, int_t end)
      : grid_(grid), start_(start), end_(end) {}

  inline Iterator begin() const { return Iterator(grid_, start_, end_); }
  inline EndIterator end() const { return EndIterator(); }

private:
  const Grid &grid_;
  int_t start_;
  int_t end_;
};

inline FaceRange exterior_faces(const Grid &grid) {
  return FaceRange(grid, grid.n_interior_edges, grid.n_edges);
}

inline FaceRange interior_faces(const Grid &grid) {
  return FaceRange(grid, 0, grid.n_interior_edges);
}

}

#endif // ZISA_FACE_RANGE_HPP
