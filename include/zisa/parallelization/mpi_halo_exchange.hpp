#ifndef ZISA_MPI_HALO_EXCHANGE_HPP_CUUIS
#define ZISA_MPI_HALO_EXCHANGE_HPP_CUUIS

#include <zisa/config.hpp>
#include <zisa/memory/array.hpp>
#include <zisa/memory/array_view.hpp>
#include <zisa/model/all_variables.hpp>
#include <zisa/parallelization/halo_exchange.hpp>
#include <zisa/parallelization/mpi.hpp>

namespace zisa {

struct HaloRemoteInfo {
  array<int_t, 1> cell_indices; ///< index local to the remote.
  int receiver_rank;
};

struct HaloLocalInfo {
  int sender_rank;

  // Locally the halo is store in the index range [i_start, i_end).
  int_t i_start;
  int_t i_end;
};

struct Halo {
  std::vector<HaloRemoteInfo> remote_info;
  std::vector<HaloLocalInfo> local_info;
};

using HaloExchangeRequest = zisa::mpi::Request;

class HaloSendPart {
private:
  using T = double;
  static constexpr int n_dims = 2;

public:
  explicit HaloSendPart(HaloRemoteInfo remote_info);

  void send(const array_const_view<T, n_dims, row_major> &out_data, int tag);

protected:
  void ensure_valid_buffer_size(const shape_t<2> &shape);
  void wait_for_send_buffer();
  void pack_buffer(const array_const_view<T, n_dims, row_major> &out_data);

private:
  HaloRemoteInfo remote_info;

  array<T, n_dims, row_major> send_buffer;
  zisa::mpi::Request send_request;
  MPI_Comm mpi_comm = MPI_COMM_WORLD;
};

class HaloReceivePart {
private:
  using T = double;
  static constexpr int n_dims = 2;

public:
  explicit HaloReceivePart(const HaloLocalInfo &local_info);

  HaloExchangeRequest receive(array_view<T, n_dims, row_major> &in_data,
                              int tag);

private:
  HaloLocalInfo local_info;
  MPI_Comm mpi_comm = MPI_COMM_WORLD;
};

class MPIHaloExchange : public HaloExchange {
private:
  using T = double;
  static constexpr int n_dims = 2;

public:
  MPIHaloExchange(std::vector<HaloRemoteInfo> remote_info,
                  const std::vector<HaloLocalInfo> &local_info);

  MPIHaloExchange(std::vector<HaloReceivePart> receive_parts,
                  std::vector<HaloSendPart> send_parts);

  void operator()(AllVariables &all_vars) override;

  void exchange(array_view<T, n_dims, row_major> data, int tag);

private:
  std::vector<HaloReceivePart> receive_parts;
  std::vector<HaloSendPart> send_parts;

  int cvars_tag = 1000;
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
std::vector<HaloRemoteInfo>
exchange_halo_info(const std::vector<HaloRemoteInfo> &remote_info,
                   const MPI_Comm &mpi_comm);

MPIHaloExchange make_mpi_halo_exchange(const Halo &halo, const MPI_Comm &comm);

}
#endif // ZISA_MPI_HALO_EXCHANGE_HPP
