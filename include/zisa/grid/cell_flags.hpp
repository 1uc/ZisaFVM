// SPDX-License-Identifier: MIT
// Copyright (c) 2021 ETH Zurich, Luc Grosheintz-Laval

#ifndef ZISA_CELL_FLAGS_HPP_UUCOO
#define ZISA_CELL_FLAGS_HPP_UUCOO

namespace zisa {

struct CellFlags {
  bool interior : 1;      ///< Cell which are 'owned'/'updated' by this grid.
  bool ghost_cell : 1;    ///< Any cell in the halo.
  bool ghost_cell_l1 : 1; ///< Has a neighbour which is an interior cell.

  CellFlags() : interior(true), ghost_cell(false), ghost_cell_l1(false) {}
};

}
#endif // ZISA_GRID_FLAGS_HPP
