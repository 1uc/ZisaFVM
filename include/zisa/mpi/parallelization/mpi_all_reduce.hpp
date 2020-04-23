#ifndef ZISA_MPI_ALL_REDUCE_HPP_CYNQP
#define ZISA_MPI_ALL_REDUCE_HPP_CYNQP

#include <zisa/config.hpp>
#include <zisa/mpi/mpi.hpp>
#include <zisa/parallelization/all_reduce.hpp>

namespace zisa {
class MPIAllReduce : public AllReduce {
public:
  MPIAllReduce(ReductionOperation op, MPI_Comm mpi_comm);

protected:
  double do_reduce(double local) const override;

private:
  MPI_Op op;
  MPI_Comm comm;
};

}
#endif // ZISA_MPI_ALL_REDUCE_HPP
