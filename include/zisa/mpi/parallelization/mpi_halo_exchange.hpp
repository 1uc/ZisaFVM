// SPDX-License-Identifier: MIT
// Copyright (c) 2021 ETH Zurich, Luc Grosheintz-Laval

#ifndef ZISA_MPI_HALO_EXCHANGE_HPP_CUUIS
#define ZISA_MPI_HALO_EXCHANGE_HPP_CUUIS

#include <map>
#include <zisa/config.hpp>
#include <zisa/memory/array.hpp>
#include <zisa/memory/array_view.hpp>
#include <zisa/model/all_variables.hpp>
#include <zisa/mpi/mpi.hpp>
#include <zisa/mpi/mpi_tag_constants.hpp>
#include <zisa/parallelization/distributed_grid.hpp>
#include <zisa/parallelization/halo_exchange.hpp>
#include <zisa/parallelization/halo_info.hpp>

namespace zisa {

using HaloExchangeRequest = zisa::mpi::Request;

class HaloSendPart {
private:
  using T = double;
  static constexpr int n_dims = 2;

public:
  explicit HaloSendPart(HaloSendInfo remote_info);

  void send(const array_const_view<T, n_dims, row_major> &out_data, int tag);

protected:
  void ensure_valid_buffer_size(const shape_t<2> &shape);
  void wait_for_send_buffer();
  void pack_buffer(const array_const_view<T, n_dims, row_major> &out_data);

private:
  HaloSendInfo send_info;

  array<T, n_dims, row_major> send_buffer;
  zisa::mpi::Request send_request;
  MPI_Comm mpi_comm = MPI_COMM_WORLD;
};

class HaloReceivePart {
private:
  using T = double;
  static constexpr int n_dims = 2;

public:
  explicit HaloReceivePart(const HaloReceiveInfo &local_info);

  HaloExchangeRequest receive(array_view<T, n_dims, row_major> &in_data,
                              int tag);

private:
  HaloReceiveInfo local_info;
  MPI_Comm mpi_comm = MPI_COMM_WORLD;
};

class MPIHaloExchange : public HaloExchange {
private:
  using T = double;
  static constexpr int n_dims = 2;

public:
  MPIHaloExchange(std::vector<HaloSendInfo> remote_info,
                  const std::vector<HaloReceiveInfo> &local_info);

  MPIHaloExchange(std::vector<HaloReceivePart> receive_parts,
                  std::vector<HaloSendPart> send_parts);

  void operator()(AllVariables &all_vars) override;

  void wait() override;

  void exchange(array_view<T, n_dims, row_major> data, int tag);

private:
  std::vector<HaloReceivePart> receive_parts;
  std::vector<HaloSendPart> send_parts;
  std::vector<HaloExchangeRequest> requests;

  int cvars_tag = ZISA_MPI_TAG_HALO_EXCHANGE_CVARS;
  int avars_tag = ZISA_MPI_TAG_HALO_EXCHANGE_AVARS;
};

/// Prepare for data exchange.
/** Problem: Rank `i` wants to send data to rank `j`, but rank `j` does
 *     not know that `i` wants to send something to it.
 *
 *  This routine accepts a vector of `(rank, size)` indicating which
 *  processors `i` wants to send data to (in bytes).
 *
 *  It returns a vector of ranks and sizes which denote which ranks want to send
 *  data to `i`.
 */
std::vector<std::pair<int, size_t>>
exchange_sizes(const std::vector<std::pair<int, size_t>> &bytes_to_send,
               const MPI_Comm &mpi_comm);

/// Exchange the remote halo information.
/** Send the information about the parts of the halo this rank needs to
 *  the remote rank which owns the data.
 *
 *  Receive the information which ranks needs what of the data owned by this
 *  rank.
 */
std::vector<HaloSendInfo>
exchange_halo_info(const std::vector<HaloRemoteInfo> &remote_info,
                   const std::map<int_t, int_t> &global2local,
                   const MPI_Comm &mpi_comm);

/// Given only the local data, create an `MPIHaloExchange`.
/** Note that this routine expects `halo.remote_info` to contain the
 *  global indices that this partition needs from the other partitions. This
 *  routine therefore internally performs a `exchange_halo_info`.
 */
std::shared_ptr<MPIHaloExchange> make_mpi_halo_exchange(
    const DistributedGrid &dgrid, const Halo &halo, const MPI_Comm &comm);

std::shared_ptr<MPIHaloExchange>
make_mpi_halo_exchange(const DistributedGrid &dgrid, const MPI_Comm &mpi_comm);

}
#endif // ZISA_MPI_HALO_EXCHANGE_HPP
