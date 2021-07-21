// SPDX-License-Identifier: MIT
// Copyright (c) 2021 ETH Zurich, Luc Grosheintz-Laval

#include <boost/program_options.hpp>

#include <zisa/grid/grid.hpp>

namespace po = boost::program_options;

int main(int argc, char *argv[]) {
  po::variables_map options;

  // generic options
  po::options_description generic("Generic options");

  // clang-format off
  generic.add_options()
      ("help,h", "produce this message")
      ("grid", po::value<std::string>(), "HDF5 grid in which a point should be located.")
      ("point", po::value<std::vector<double>>(), "Coordinates of the point to locate.")
      ;

  po::positional_options_description positional;
  positional.add("point", -1);

  po::store(
      po::command_line_parser(argc, argv).options(generic)
                                         .positional(positional)
                                         .run(),
      options);
  po::notify(options);

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

  if (options.count("point") == 0) {
    std::cout << "Missing coordinate of the point.\n";
    std::exit(EXIT_FAILURE);
  }

  const auto &x_ = options["point"].as<std::vector<double>>();
  auto x = (x_.size() == 2 ? zisa::XYZ{x_[0], x_[1]}
                           : zisa::XYZ{x_[0], x_[1], x_[2]});

  auto grid_file = options["grid"].as<std::string>();
  auto grid = zisa::load_grid(grid_file);

  auto i = zisa::locate(*grid, x);

  if (i == std::nullopt) {
    std::cerr << "Point not found.\n";
    std::exit(EXIT_FAILURE);
  }

  std::cout << i.value() << "\n";
  std::exit(EXIT_SUCCESS);
}
