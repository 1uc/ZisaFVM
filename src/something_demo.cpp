#include <iostream>
#include <string>

#include <zisa/model/helmholtz_eos.hpp>
#include <zisa/reconstruction/cweno_ao.hpp>

namespace zisa {

void test_helmholtz_eos() {
  auto table_path = std::string("data/stellar_convection/helm_table.dat");
  initialize_helmholtz_eos(table_path);

  //  auto mass_mixing_ratio = std::vector<double>{
  //      0.0013448965927533212,
  //      9.39876715020613e-21,
  //      6.132204407482413e-15,
  //      0.3949170331615017,
  //      0.5288770106729471,
  //      0.0021922858917659254,
  //  };
  //  auto mass_number = std::vector<double>{12, 1, 4, 20, 16, 28};
  //  auto charge_number = std::vector<double>{6, 1, 2, 10, 8, 14};

  //  double h = 7.853e+16;
  //  double s = 4.566e+08;
  //  double rho = 5.122e+03;
  //  double T = 3.638e+08;
  //  double a_bar = 1.484027e+01;
  //  double z_bar = 7.278513e+00;

  double e = 5.033330e+14;
  double rho = 1.251521e+03;
  double a_bar = 1.484026e+01;
  double z_bar = 7.278644e+00;

  auto eos = HelmholtzEOS(a_bar, z_bar);

  //  double e = 5.335833e+16;
  //  double rho = 2.464040e+03;

  auto full_xvars = eos.full_extra_variables(RhoE{rho, rho * e});
  //  auto full_xvars
  //      = eos.full_extra_variables(EnthalpyEntropy{h, s}, RhoT{rho, T});
  PRINT(full_xvars.p);
  PRINT(full_xvars.E);
  PRINT(full_xvars.T);
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
