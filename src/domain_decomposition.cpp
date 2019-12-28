#include <iostream>
#include <string>
#include <boost/program_options.hpp>

#include <metis.h>

#include <zisa/grid/grid.hpp>
#include <zisa/io/format_as_list.hpp>

namespace po = boost::program_options;

zisa::array<idx_t, 1> compute_partitions(const zisa::Grid &grid, int n_parts) {
  auto n_cells = grid.n_cells;
  auto n_vertices = grid.n_vertices;
  auto max_neighbours = grid.max_neighbours;
  const auto &vertex_indices = grid.vertex_indices;

  auto ne = idx_t(n_cells);
  auto nn = idx_t(n_vertices);

  zisa::array<idx_t, 1> eptr(n_cells+1);
  zisa::array<idx_t, 1> eind({n_cells * max_neighbours});

  for(zisa::int_t i = 0; i < n_cells; ++i) {
    eptr[i] = 3*idx_t(i);

    for(zisa::int_t k = 0; k < max_neighbours; ++k) {
      eind[eptr[i] + k] = vertex_indices(i, k);
    }
  }
  eptr[n_cells] = 3*idx_t(n_cells);

  auto ncommon = idx_t(max_neighbours);
  auto nparts = idx_t(n_parts);

  idx_t objval;
  zisa::array<idx_t, 1> epart({n_cells});
  zisa::array<idx_t, 1> npart({n_vertices});

  METIS_PartMeshDual(&ne, &nn, eptr.raw(), eind.raw(), nullptr, nullptr,
                     &ncommon, &nparts, nullptr, nullptr, &objval, epart.raw(), npart.raw());

  return epart;
}

int main(int argc, char *argv[]) {
  po::variables_map options;

  // generic options
  po::options_description generic("Generic options");

  // clang-format off
  generic.add_options()
      ("help,h", "produce this message")
      ("grid", po::value<std::string>(), "GMSH grid to partition.")
      ("partitions,n", po::value<int>(), "Number of partitions.")
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

  if (options.count("partitions") == 0) {
    std::cout << "Missing argument `-n N`.\n";
    std::exit(EXIT_FAILURE);
  }

  auto n_parts = options["partitions"].as<int>();
  auto gmsh_file = options["grid"].as<std::string>();
  auto grid = zisa::load_gmsh(gmsh_file);

  std::cout << zisa::format_as_list(compute_partitions(*grid, n_parts)) << "\n";
}
