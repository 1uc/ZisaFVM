#include <iostream>
#include <string>

#include <boost/program_options.hpp>
namespace po = boost::program_options;

#include <filesystem>
#include <zisa/io/hdf5_serial_writer.hpp>
#include <zisa/math/cartesian.hpp>
#include <zisa/math/permutation.hpp>
#include <zisa/math/space_filling_curve.hpp>
#include <zisa/memory/array.hpp>
#include <zisa/memory/array_view.hpp>

namespace zisa {
void renumber_grid(const std::string &grid_file) {
  auto [vertices, vertex_indices, n_dims] = [&grid_file]() {
    auto reader = HDF5SerialReader(grid_file);

    auto vertices = array<XYZ, 1>::load(reader, "vertices");
    auto vertex_indices = array<int_t, 2>::load(reader, "vertex_indices");
    auto n_dims = reader.read_scalar<int_t>("n_dims");

    return std::tuple{
        std::move(vertices), std::move(vertex_indices), std::move(n_dims)};
  }();

  auto n_cells = vertex_indices.shape(0);
  auto max_vertices = vertex_indices.shape(1);
  auto cell_centers = array<XYZ, 1>(n_cells);

  double x_min = std::numeric_limits<double>::max();
  double x_max = std::numeric_limits<double>::min();

  double y_min = std::numeric_limits<double>::max();
  double y_max = std::numeric_limits<double>::min();

  double z_min = std::numeric_limits<double>::max();
  double z_max = std::numeric_limits<double>::min();

  for (int_t i = 0; i < n_cells; ++i) {
    cell_centers[i] = vertices[vertex_indices(i, 0)];
    for (int_t k = 0; k < max_vertices; ++k) {
      cell_centers[i] += vertices[vertex_indices(i, k)];
    }

    cell_centers[i] /= max_vertices;

    x_min = zisa::min(cell_centers[i][0], x_min);
    x_max = zisa::max(cell_centers[i][0], x_max);

    y_min = zisa::min(cell_centers[i][1], y_min);
    y_max = zisa::max(cell_centers[i][1], y_max);

    z_min = zisa::min(cell_centers[i][2], z_min);
    z_max = zisa::max(cell_centers[i][2], z_max);
  }

  for (int_t i = 0; i < n_cells; ++i) {
    auto [x, y, z] = cell_centers[i];

    cell_centers[i][0] = (x - x_min) / (x_max - x_min + 1e-10 * x_max);
    cell_centers[i][1] = (y - y_min) / (y_max - y_min + 1e-10 * y_max);

    if (n_dims == 3) {
      cell_centers[i][2] = (z - z_min) / (z_max - z_min + 1e-10 * z_max);
    }
  }

  auto sfc_indices = array<int_t, 1>(n_cells);
  for (int_t i = 0; i < n_cells; ++i) {
    if (n_dims == 2) {
      auto [x, y, _] = cell_centers[i];
      sfc_indices[i]
          = integer_cast<int_t>(hilbert_index<64 / 2>(x, y).to_ullong());
    }

    if (n_dims == 3) {
      auto [x, y, z] = cell_centers[i];
      sfc_indices[i]
          = integer_cast<int_t>(hilbert_index<64 / 3>(x, y, z).to_ullong());
    }
  }

  auto sigma = array<int_t, 1>(n_cells);
  for (int_t i = 0; i < n_cells; ++i) {
    sigma[i] = i;
  }
  std::sort(sigma.begin(), sigma.end(), [&sfc_indices](int_t i, int_t j) {
    return sfc_indices[i] < sfc_indices[j];
  });

  apply_permutation(array_view(vertex_indices), factor_permutation(sigma));

  {
    auto writer = HDF5SerialWriter(grid_file + "_");
    save(writer, vertices, "vertices");
    save(writer, vertex_indices, "vertex_indices");
    writer.write_scalar(n_dims, "n_dims");
  }
  std::filesystem::rename(grid_file + "_", grid_file);
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
      ("grid", po::value<std::string>(), "Name of the grid file, will be overwritten.")
      ;
  // clang-format on

  // first parse cmdline and check what config file to use
  po::store(po::parse_command_line(argc, argv, generic), options);

  if (options.count("help") != 0) {
    std::cout << generic << "\n";
    std::exit(EXIT_SUCCESS);
  }

  if (options.count("grid") == 0) {
    std::cout << "Missing argument `--grid GRID`.\n";
    std::exit(EXIT_FAILURE);
  }

  auto grid_file = options["grid"].as<std::string>();

  zisa::renumber_grid(grid_file);

  MPI_Finalize();
}
