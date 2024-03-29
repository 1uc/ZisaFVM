// SPDX-License-Identifier: MIT
// Copyright (c) 2021 ETH Zurich, Luc Grosheintz-Laval

#include <zisa/mpi/io/mpi_progress_bar.hpp>

namespace zisa {

MPIProgressBar::MPIProgressBar(std::shared_ptr<ProgressBar> serial_bar,
                               MPI_Comm mpi_comm)
    : serial_bar(std::move(serial_bar)), rank(zisa::mpi::rank(mpi_comm)) {}

void MPIProgressBar::write_progress(std::ostream &out, const std::string &msg) {
  if (rank == 0) {
    serial_bar->write_progress(out, msg);
  }
}

void MPIProgressBar::activate() { serial_bar->activate(); }
void MPIProgressBar::deactivate() { serial_bar->deactivate(); }
void MPIProgressBar::reset() { serial_bar->reset(); }

}
