#define CATCH_CONFIG_RUNNER
#include "zisa/testing/testing_framework.hpp"

#if ZISA_HAS_MPI == 1
#include <zisa/mpi/mpi.hpp>
#endif

int main( int argc, char* argv[] ) {
#if ZISA_HAS_MPI == 1
  int requested = MPI_THREAD_MULTIPLE;
  int provided = -1;

  MPI_Init_thread(&argc, &argv, requested, &provided);
  LOG_ERR_IF(requested > provided,
             "MPI does not support the requested level of multi-threading.");
#endif

  int result = Catch::Session().run( argc, argv );

#if ZISA_HAS_MPI == 1
  MPI_Finalize();
#endif

  return result;
}