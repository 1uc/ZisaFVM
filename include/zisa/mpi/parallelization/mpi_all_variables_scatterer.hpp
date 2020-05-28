#ifndef ZISA_MPI_ALL_VARIABLES_SCATTERER_HPP_XBTQI
#define ZISA_MPI_ALL_VARIABLES_SCATTERER_HPP_XBTQI

#include <zisa/mpi/mpi.hpp>
#include <zisa/parallelization/all_variables_scatterer.hpp>
#include <zisa/memory/array.hpp>

namespace {

std::unique_ptr<AllVariablesScatterer> make_mpi_all_variables_scatterer(
    int base_tag, const MPI_Comm &comm, const array<int_t, 1> &boundaries);

}

#endif // ZISA_MPI_ALL_VARIABLES_SCATTERER_HPP
