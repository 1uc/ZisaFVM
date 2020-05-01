#ifndef ZISA_MPI_ALL_VARIABLES_GATHERER_HPP
#define ZISA_MPI_ALL_VARIABLES_GATHERER_HPP

#include <zisa/mpi/mpi.hpp>
#include <zisa/mpi/parallelization/mpi_single_node_array_gatherer.hpp>
#include <zisa/parallelization/all_variables_gatherer.hpp>

namespace zisa {

std::unique_ptr<AllVariablesGatherer> make_mpi_all_variables_gatherer(
    int base_tag, const MPI_Comm &comm, const array<int_t, 1> &boundaries);

}

#endif // ZISA_MPI_ALL_VARIABLES_GATHERER_HPP