// SPDX-License-Identifier: MIT
// Copyright (c) 2021 ETH Zurich, Luc Grosheintz-Laval

#include <zisa/loops/range.hpp>

#include <zisa/model/grid_variables.hpp>

namespace zisa {

Range<PlainIndexRange, DereferenceConstGridVariables>
cell_const_range(const GridVariables &grid_vars) {
  return Range(PlainIndexRange(0, grid_vars.shape(0)),
               DereferenceConstGridVariables(grid_vars));
}

}