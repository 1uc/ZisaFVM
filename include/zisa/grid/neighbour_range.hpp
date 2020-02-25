#ifndef ZISA_NEIGHBOUR_RANGE_HPP_IUCKK
#define ZISA_NEIGHBOUR_RANGE_HPP_IUCKK

#include <zisa/config.hpp>

#include <zisa/grid/grid.hpp>
#include <zisa/loops/range.hpp>

namespace zisa {

inline PlainIndexRange neighbour_index_range(const Grid &grid) {
  return PlainIndexRange(0, grid.max_neighbours);
}

}

#endif // ZISA_NEIGHBOUR_RANGE_HPP
