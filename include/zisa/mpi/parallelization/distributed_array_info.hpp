// SPDX-License-Identifier: MIT
// Copyright (c) 2021 ETH Zurich, Luc Grosheintz-Laval

#ifndef ZISA_DISTRIBUTED_ARRAY_INFO_HPP_ZXXKX
#define ZISA_DISTRIBUTED_ARRAY_INFO_HPP_ZXXKX

#include <zisa/config.hpp>

#include <zisa/memory/array.hpp>

namespace zisa {

struct DistributedArrayInfo {
  explicit DistributedArrayInfo(array<int_t, 1> partition);

  // Rank `p` has all data with indices [partition[p], partition[p+1]).
  array<int_t, 1> partition;
};

std::shared_ptr<DistributedArrayInfo>
make_distributed_array_info(array<int_t, 1> partition);

}

#endif // ZISA_DISTRIBUTED_ARRAY_INFO_HPP
