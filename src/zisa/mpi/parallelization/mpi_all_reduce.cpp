#include <zisa/mpi/parallelization/mpi_all_reduce.hpp>

namespace zisa {

MPIAllReduce::MPIAllReduce(ReductionOperation op, MPI_Comm mpi_comm)
    : op(MPI_MIN), comm(mpi_comm) {

  LOG_ERR_IF(op != ReductionOperation::min, "Implement the other cases.");
}

double MPIAllReduce::do_reduce(double local) const {
  double global = 0.0;

  auto code = MPI_Allreduce(&local, &global, 1, MPI_DOUBLE, op, comm);
  LOG_ERR_IF(code != MPI_SUCCESS,
             string_format("MPI_Allreduce failed. [%d]", code));

  return global;
}

}