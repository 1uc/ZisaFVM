#import <zisa/mpi/math/distributed_reference_solution.hpp>

namespace zisa {

void test_distributed_reference() {

  //  auto dref = DistributedReferenceSolution(
  //      serialize, large_grid, interpolation, n_vars);
}

}

int main(int argc, char *argv[]) {
  MPI_Init(&argc, &argv);

  zisa::test_distributed_reference();

  MPI_Finalize();
}
