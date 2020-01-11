#ifndef ZISA_HALO_EXCHANGE_HPP_IDOIW
#define ZISA_HALO_EXCHANGE_HPP_IDOIW

#include <zisa/config.hpp>
#include <zisa/memory/array.hpp>
#include <zisa/model/all_variables.hpp>
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

using HaloExchangeRequest = zisa::mpi::Request;

class HaloSendPart {
private:
  using T = double;
  static constexpr int n_dims = 2;

public:
  explicit HaloSendPart(HaloRemoteInfo remote_info)
      : remote_info(std::move(remote_info)), send_request(nullptr) {}

  void send(const array_const_view<T, n_dims, row_major> &out_data, int tag) {
    wait_for_send_buffer();
    ensure_valid_buffer_size(out_data.shape());
    pack_buffer(out_data);

    auto receiver = remote_info.receiver_rank;
    auto const_view = array_const_view(send_buffer);
    send_request = zisa::mpi::isend(const_view, receiver, tag, mpi_comm);
  }

protected:
  void ensure_valid_buffer_size(const shape_t<2> &shape) {
    int_t n_cells_buffer = remote_info.cell_indices.shape(0);

    if (send_buffer.shape(0) != n_cells_buffer
        || send_buffer.shape(1) != shape(1)) {
      send_buffer = array<T, n_dims, row_major>({n_cells_buffer, shape(1)});
    }
  }

  void wait_for_send_buffer() {
    send_request.wait();
  }

  void pack_buffer(const array_const_view<T, n_dims, row_major> &out_data) {
    const auto &cell_indices = remote_info.cell_indices;

    for (int_t i = 0; i < cell_indices.size(); ++i) {
      for (int_t k = 0; k < out_data.shape(1); ++k) {
        send_buffer(i, k) = out_data(cell_indices(i), k);
      }
    }
  }

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
  explicit HaloReceivePart(const HaloLocalInfo &local_info)
      : local_info(local_info) {}

  HaloExchangeRequest receive(array_view<T, n_dims, row_major> &in_data,
                              int tag) {
    auto sub_array = slice(in_data, local_info.i_start, local_info.i_end);
    return zisa::mpi::irecv(sub_array, local_info.sender_rank, tag, mpi_comm);
  }

private:
  HaloLocalInfo local_info;
  MPI_Comm mpi_comm = MPI_COMM_WORLD;
};

class HaloExchange {
public:
  virtual void operator()(AllVariables &all_vars) = 0;
};

class MPIHaloExchange : public HaloExchange {
private:
  using T = double;
  static constexpr int n_dims = 2;

public:
  MPIHaloExchange(std::vector<HaloRemoteInfo> remote_info,
                  std::vector<HaloLocalInfo> local_info) {
    send_parts.reserve(remote_info.size());
    for (auto &r : remote_info) {
      send_parts.emplace_back(std::move(r));
    }

    receive_parts.reserve(local_info.size());
    for (auto &l : local_info) {
      receive_parts.emplace_back(std::move(l));
    }
  }

  MPIHaloExchange(std::vector<HaloReceivePart> receive_parts,
                  std::vector<HaloSendPart> send_parts)
      : receive_parts(std::move(receive_parts)),
        send_parts(std::move(send_parts)) {}

  void operator()(AllVariables &all_vars) override {
    LOG_ERR_IF(
        all_vars.avars.shape(1) != 0,
        "Halo exchange of advected variables need to be implemented first.");

    exchange(all_vars.cvars, cvars_tag);
  }

  void exchange(array_view<T, n_dims, row_major> data, int tag) {
    std::vector<HaloExchangeRequest> requests;
    requests.reserve(receive_parts.size());

    for (auto &p : receive_parts) {
      requests.push_back(p.receive(data, tag));
    }

    for (auto &p : send_parts) {
      p.send(data, tag);
    }

    for (const auto &r : requests) {
      r.wait();
    }
  }

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
               const MPI_Comm &mpi_comm) {
  auto mpi_rank = zisa::mpi::rank(mpi_comm);
  auto mpi_ranks = zisa::mpi::size(mpi_comm);
  std::vector<size_t> send_buffer(mpi_ranks, size_t(0));

  for (const auto &[rank, size] : bytes_to_send) {
    send_buffer[rank] = size;
  }

  std::vector<size_t> recv_buffer(mpi_ranks);
  for (int i = 0; i < mpi_ranks; ++i) {
    if (i != mpi_rank) {
      auto status = MPI_Scatter(nullptr,
                                0,
                                MPI_BYTE,
                                (void *)&recv_buffer[i],
                                sizeof(size_t),
                                MPI_BYTE,
                                i,
                                mpi_comm);
      LOG_ERR_IF(status != MPI_SUCCESS,
                 string_format(
                     "MPI_Scatter failed (rank = %d). [%d]", mpi_rank, status));
    } else {
      auto status = MPI_Scatter(send_buffer.data(),
                                sizeof(size_t),
                                MPI_BYTE,
                                (void *)&recv_buffer[i],
                                sizeof(size_t),
                                MPI_BYTE,
                                i,
                                mpi_comm);
      LOG_ERR_IF(status != MPI_SUCCESS,
                 string_format(
                     "MPI_Scatter failed (root = %d). [%d]", mpi_rank, status));
    }
  }

  std::vector<std::pair<int, size_t>> bytes_to_receive;
  for (int i = 0; i < mpi_ranks; ++i) {
    if (recv_buffer[i] != 0) {
      bytes_to_receive.emplace_back(i, recv_buffer[i]);
    }
  }

  return bytes_to_receive;
}

/// Exchange the remote halo information.
/** Send the information about the parts of the halo this rank needs to
 *  the remote rank which owns the data.
 *
 *  Receive the information which ranks needs what of the data owned by this
 *  rank.
 */
std::vector<HaloRemoteInfo>
exchange_halo_info(const std::vector<HaloRemoteInfo> &remote_info,
                   const MPI_Comm &mpi_comm) {
  int mpi_rank = zisa::mpi::rank(mpi_comm);
  int mpi_ranks = zisa::mpi::size(mpi_comm);

  std::vector<std::pair<int, size_t>> bytes_to_send;
  bytes_to_send.reserve(remote_info.size());

  for (const auto &r : remote_info) {
    auto rank = r.receiver_rank;
    auto size = r.cell_indices.size() * sizeof(int_t);
    bytes_to_send.emplace_back(rank, size);
  }

  auto bytes_to_receive = exchange_sizes(bytes_to_send, mpi_comm);

  int xfer_tag = 200;

  // Post all receives
  std::vector<HaloRemoteInfo> received_remote_info;
  std::vector<MPI_Request> recv_requests;
  for (auto [rank, size] : bytes_to_receive) {
    recv_requests.emplace_back();

    int_t n_cells_halo = size / sizeof(int_t);
    received_remote_info.emplace_back(
        HaloRemoteInfo{array<int_t, 1>(n_cells_halo), rank});
    const auto &cell_indices = received_remote_info.back().cell_indices;

    MPI_Irecv((void *)cell_indices.raw(),
              size,
              MPI_BYTE,
              rank,
              xfer_tag,
              mpi_comm,
              &recv_requests.back());
  }

  std::vector<MPI_Request> send_requests;
  for (const auto &r : remote_info) {
    auto ptr = (void *)r.cell_indices.raw();
    auto size = r.cell_indices.size() * sizeof(int_t);

    send_requests.emplace_back();
    MPI_Isend(ptr,
              size,
              MPI_BYTE,
              r.receiver_rank,
              xfer_tag,
              mpi_comm,
              &send_requests.back());
  }

  MPI_Waitall(send_requests.size(), send_requests.data(), MPI_STATUSES_IGNORE);
  MPI_Waitall(recv_requests.size(), recv_requests.data(), MPI_STATUSES_IGNORE);

  return received_remote_info;
}

}

#endif // ZISA_HALO_EXCHANGE_HPP
