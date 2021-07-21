// SPDX-License-Identifier: MIT
// Copyright (c) 2021 ETH Zurich, Luc Grosheintz-Laval

#ifndef ZISA_HALO_INFO_HPP_CNUWI
#define ZISA_HALO_INFO_HPP_CNUWI

#include <zisa/config.hpp>

#include <vector>
#include <zisa/memory/array.hpp>

namespace zisa {

/// This describes what this PE needs from `remote_rank`.
struct HaloRemoteInfo {
  int remote_rank;
  array<int_t, 1> cell_indices; ///< these are global indices.

  HaloRemoteInfo(int remote_rank, array<int_t, 1> cell_indices)
      : remote_rank(remote_rank), cell_indices(std::move(cell_indices)) {}
};

/// This describes what `receiver_rank` needs from us.
struct HaloSendInfo {
  int receiver_rank;
  array<int_t, 1> cell_indices; ///< index local to this PE.

  HaloSendInfo(int receiver_rank, array<int_t, 1> cell_indices)
      : receiver_rank(receiver_rank), cell_indices(std::move(cell_indices)) {}
};

/// This describes the halo for `sender_rank` on the local PE.
struct HaloReceiveInfo {
  int sender_rank;

  // Locally the halo is stored in the index range [i_start, i_end).
  int_t i_start;
  int_t i_end;

  HaloReceiveInfo(int sender_rank, int_t i_start, int_t i_end)
      : sender_rank(sender_rank), i_start(i_start), i_end(i_end) {}
};

struct Halo {
  std::vector<HaloRemoteInfo> remote_info;
  std::vector<HaloReceiveInfo> local_info;
};

}

#endif // ZISA_HALO_INFO_HPP
