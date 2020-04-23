#if ZISA_HAS_MPI == 1
#include <iostream>
#include <string>
#include <tuple>

#include <zisa/mpi/io/hdf5_unstructured_writer.hpp>
#include <zisa/mpi/mpi.hpp>

#include <boost/program_options.hpp>
namespace po = boost::program_options;

namespace zisa {
void test_parallel_hdf5(const std::string &output_file) {
  LOG_ERR_IF(zisa::mpi::size() != 2, "Must be run with 2 MPI ranks.");

  auto mpi_rank = zisa::mpi::rank();

  auto ids = (mpi_rank == 0 ? std::vector<hsize_t>{0, 2, 3, 9}
                            : std::vector<hsize_t>{1, 4, 5, 6, 7, 8});

  auto upsilon = array<double, 1>(ids.size());
  for (auto i : index_range(ids.size())) {
    upsilon(i) = ids[i];
  }

  auto tau = array<double, 2>(shape_t<2>{ids.size() + 4, 3});
  for (auto i : index_range(ids.size())) {
    tau(i, 0) = 100 * ids[i];
    tau(i, 1) = 100 * ids[i] + 1;
    tau(i, 2) = 100 * ids[i] + 2;
  }

  auto file_dims = make_hdf5_unstructured_file_dimensions(ids, MPI_COMM_WORLD);

  auto writer = HDF5UnstructuredWriter(output_file, file_dims);
  save(writer, upsilon, "upsilon");
  save(writer, tau, "tau");
}
}

int main(int argc, char *argv[]) {
  MPI_Init(&argc, &argv);

  po::variables_map options;

  // generic options
  po::options_description generic("Generic options");

  // clang-format off
  generic.add_options()
      ("help,h", "produce this message")
      ("output,o", po::value<std::string>(), "Name of the output file.")
      ;
  // clang-format on

  // first parse cmdline and check what config file to use
  po::store(po::parse_command_line(argc, argv, generic), options);

  if (options.count("help") != 0) {
    std::cout << generic << "\n";
    std::exit(EXIT_SUCCESS);
  }

  if (options.count("output") == 0) {
    std::cout << "Missing argument `--output OUTPUT`.\n";
    std::exit(EXIT_FAILURE);
  }

  auto output_file = options["output"].as<std::string>();

  zisa::test_parallel_hdf5(output_file);

  MPI_Finalize();
  return EXIT_SUCCESS;
}

#else
#include <zisa/config.hpp>
int main(int argc, char *argv[]) { LOG_ERR("Needs to be compiles with MPI"); }
#endif
