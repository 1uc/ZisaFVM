#include <iostream>
#include <random>
#include <string>
#include <tuple>

#include <Eigen/Dense>
#include <zisa/model/helmholtz_eos.hpp>
#include <zisa/reconstruction/cweno_ao.hpp>

namespace zisa {

void test_helmholtz_eos() {
  auto eos = HelmholtzEOS("foo");
  auto cweno = CWENO_AO{};
}

}

int main(int argc, char *argv[]) {
#if ZISA_HAS_MPI == 1
  int available_thread_level = -1;
  int requested_thread_level = MPI_THREAD_MULTIPLE;
  MPI_Init_thread(
      &argc, &argv, requested_thread_level, &available_thread_level);
  LOG_ERR_IF(available_thread_level != requested_thread_level,
             "Can't init MPI.");
#endif

  zisa::test_helmholtz_eos();

#if ZISA_HAS_MPI == 1
  MPI_Finalize();
#endif
  return EXIT_SUCCESS;
}
