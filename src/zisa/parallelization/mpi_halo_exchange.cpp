#include <zisa/parallelization/mpi_halo_exchange.hpp>

#include <zisa/io/format_as_list.hpp>
#include <zisa/parallelization/distributed_grid.hpp>

namespace zisa {

std::vector<std::pair<int, size_t>>
exchange_sizes(const std::vector<std::pair<int, size_t>> &bytes_to_send,
               const MPI_Comm &mpi_comm) {
  auto mpi_rank = zisa::mpi::rank(mpi_comm);
  auto mpi_ranks = zisa::mpi::size(mpi_comm);

  auto n_ranks = integer_cast<size_t>(mpi_ranks);

  std::vector<size_t> send_buffer(n_ranks, size_t(0));

  for (const auto &[rank, size] : bytes_to_send) {
    send_buffer[integer_cast<size_t>(rank)] = size;
  }

  std::vector<size_t> recv_buffer(n_ranks);
  for (size_t i = 0; i < n_ranks; ++i) {
    auto status = MPI_Scatter(send_buffer.data(),
                              sizeof(size_t),
                              MPI_BYTE,
                              (void *)&recv_buffer[i],
                              sizeof(size_t),
                              MPI_BYTE,
                              integer_cast<int>(i),
                              mpi_comm);
    LOG_ERR_IF(status != MPI_SUCCESS,
               string_format(
                   "MPI_Scatter failed (root = %d). [%d]", mpi_rank, status));
  }

  std::vector<std::pair<int, size_t>> bytes_to_receive;
  for (size_t i = 0; i < n_ranks; ++i) {
    if (recv_buffer[i] != 0) {
      bytes_to_receive.emplace_back(i, recv_buffer[i]);
    }
  }

  return bytes_to_receive;
}

std::vector<HaloSendInfo>
exchange_halo_info(const std::vector<HaloRemoteInfo> &remote_info,
                   const std::map<int_t, int_t> &g2l,
                   const MPI_Comm &mpi_comm) {

  std::vector<std::pair<int, size_t>> bytes_to_send;
  bytes_to_send.reserve(remote_info.size());

  for (const auto &r : remote_info) {
    auto rank = r.remote_rank;
    auto size = r.cell_indices.size() * sizeof(int_t);
    bytes_to_send.emplace_back(rank, size);
  }

  auto bytes_to_receive = exchange_sizes(bytes_to_send, mpi_comm);

  int xfer_tag = ZISA_MPI_TAG_EXCHANGE_HALO_INFO_XFER;

  // Post all receives
  std::vector<HaloSendInfo> send_info;
  std::vector<zisa::mpi::Request> recv_requests;
  recv_requests.reserve(bytes_to_receive.size());
  for (auto [rank, size] : bytes_to_receive) {
    int_t n_cells_halo = size / sizeof(int_t);
    send_info.emplace_back(HaloSendInfo{rank, array<int_t, 1>(n_cells_halo)});
    auto &cell_indices = send_info.back().cell_indices;

    recv_requests.push_back(
        zisa::mpi::irecv(array_view(cell_indices), rank, xfer_tag, mpi_comm));
  }

  std::vector<zisa::mpi::Request> send_requests;
  send_requests.reserve(send_requests.size());
  for (const auto &r : remote_info) {
    send_requests.push_back(zisa::mpi::isend(
        array_const_view(r.cell_indices), r.remote_rank, xfer_tag, mpi_comm));
  }

  zisa::mpi::wait_all(send_requests);
  zisa::mpi::wait_all(recv_requests);

  // Convert from global to local indices.
  for (auto &si : send_info) {
    for (auto &i : si.cell_indices) {
      i = g2l.at(i);
    }
  }

  return send_info;
}

HaloSendPart::HaloSendPart(HaloSendInfo remote_info)
    : send_info(std::move(remote_info)), send_request(nullptr) {}

void HaloSendPart::send(const array_const_view<T, n_dims, row_major> &out_data,
                        int tag) {
  wait_for_send_buffer();
  ensure_valid_buffer_size(out_data.shape());
  pack_buffer(out_data);

  auto receiver = send_info.receiver_rank;
  auto const_view = array_const_view(send_buffer);

  send_request = zisa::mpi::isend(const_view, receiver, tag, mpi_comm);
}

void HaloSendPart::ensure_valid_buffer_size(const shape_t<2> &shape) {
  int_t n_cells_buffer = send_info.cell_indices.shape(0);

  if (send_buffer.shape(0) != n_cells_buffer
      || send_buffer.shape(1) != shape(1)) {
    send_buffer = array<T, n_dims, row_major>({n_cells_buffer, shape(1)});
  }
}

void HaloSendPart::wait_for_send_buffer() { send_request.wait(); }

void HaloSendPart::pack_buffer(
    const array_const_view<T, n_dims, row_major> &out_data) {

  const auto &cell_indices = send_info.cell_indices;

  for (int_t i = 0; i < cell_indices.size(); ++i) {
    for (int_t k = 0; k < out_data.shape(1); ++k) {
      send_buffer(i, k) = out_data(cell_indices(i), k);
    }
  }
}

HaloReceivePart::HaloReceivePart(const HaloReceiveInfo &local_info)
    : local_info(local_info) {}

HaloExchangeRequest
HaloReceivePart::receive(array_view<T, n_dims, row_major> &in_data, int tag) {
  auto sub_array = slice(in_data, local_info.i_start, local_info.i_end);
  return zisa::mpi::irecv(sub_array, local_info.sender_rank, tag, mpi_comm);
}

MPIHaloExchange::MPIHaloExchange(
    std::vector<HaloSendInfo> remote_info,
    const std::vector<HaloReceiveInfo> &local_info) {

  send_parts.reserve(remote_info.size());
  for (auto &r : remote_info) {
    send_parts.emplace_back(std::move(r));
  }

  receive_parts.reserve(local_info.size());
  for (const auto &l : local_info) {
    receive_parts.emplace_back(l);
  }
}

MPIHaloExchange::MPIHaloExchange(std::vector<HaloReceivePart> receive_parts,
                                 std::vector<HaloSendPart> send_parts)
    : receive_parts(std::move(receive_parts)),
      send_parts(std::move(send_parts)) {}

void MPIHaloExchange::operator()(AllVariables &all_vars) {
  LOG_ERR_IF(
      all_vars.avars.shape(1) != 0,
      "Halo exchange of advected variables need to be implemented first.");

  exchange(all_vars.cvars, cvars_tag);
}

void MPIHaloExchange::exchange(array_view<T, n_dims, row_major> data, int tag) {
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

MPIHaloExchange make_mpi_halo_exchange(const DistributedGrid &dgrid,
                                       const Halo &halo,
                                       const MPI_Comm &comm) {
  auto g2l = make_global2local(dgrid.global_cell_indices);
  auto remote_info = exchange_halo_info(halo.remote_info, g2l, comm);
  auto &local_info = halo.local_info;

  return zisa::MPIHaloExchange(std::move(remote_info), local_info);
}

}
