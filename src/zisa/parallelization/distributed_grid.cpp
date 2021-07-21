// SPDX-License-Identifier: MIT
// Copyright (c) 2021 ETH Zurich, Luc Grosheintz-Laval

#include <zisa/parallelization/distributed_grid.hpp>

#include <zisa/io/hdf5_serial_writer.hpp>

namespace zisa {

DistributedGrid load_distributed_grid(const std::string &filename) {
  auto reader = HDF5SerialReader(filename);

  auto gci = array<int_t, 1>::load(reader, "global_cell_indices");
  auto p = array<int_t, 1>::load(reader, "partition");

  return DistributedGrid{std::move(gci), std::move(p)};
}

std::map<int_t, int_t> make_global2local(const array<int_t, 1> &l2g) {

  auto g2l = std::map<int_t, int_t>();

  auto n_cells = l2g.size();
  for (int_t i = 0; i < n_cells; ++i) {
    g2l[l2g[i]] = i;
  }

  g2l[int_t(-1)] = int_t(-1);

  return g2l;
}
}