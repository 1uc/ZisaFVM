#ifndef ZISA_DISTRIBUTED_GRID_INFO_HPP_ZNASJ
#define ZISA_DISTRIBUTED_GRID_INFO_HPP_ZNASJ

#include <zisa/config.hpp>

#include <map>
#include <zisa/memory/array.hpp>

namespace zisa {

struct DistributedGrid {
  array<int_t, 1> global_cell_indices;
  array<int_t, 1> partition;
};

std::map<int_t, int_t> make_global2local(const array<int_t, 1> &l2g);

DistributedGrid load_distributed_grid(const std::string &filename);

}
#endif // ZISA_DISTRIBUTED_GRID_HPP
