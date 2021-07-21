// SPDX-License-Identifier: MIT
// Copyright (c) 2021 ETH Zurich, Luc Grosheintz-Laval

#ifndef ZISA_ARRAY_CELL_FLAGS_HPP
#define ZISA_ARRAY_CELL_FLAGS_HPP

#include <zisa/grid/cell_flags.hpp>
#include <zisa/memory/array.hpp>

namespace zisa {

struct cell_flags_dispatch_tag {};

template <>
struct array_save_traits<CellFlags> {
  using dispatch_tag = cell_flags_dispatch_tag;
};

template <int n_dims>
void save(HierarchicalWriter &writer,
          const array_const_view<CellFlags, n_dims, row_major> &arr,
          const std::string &tag,
          cell_flags_dispatch_tag) {

  using scalar_type = typename array_save_traits<bool>::scalar_type;
  auto int_arr = array<scalar_type, n_dims>(arr.shape());

  writer.open_group(tag);
  std::transform(
      arr.cbegin(),
      arr.cend(),
      int_arr.begin(),
      [](const CellFlags &cell_flags) { return cell_flags.ghost_cell; });
  save(writer, int_arr, "ghost_cell");

  writer.close_group();
}

}
#endif // ZISA_ARRAY_CELL_FLAGS_HPP
