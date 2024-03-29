// SPDX-License-Identifier: MIT
// Copyright (c) 2021 ETH Zurich, Luc Grosheintz-Laval

/* List of grids that can be used for convergence studies.
 */

#ifndef GRIDS_H_EUWJG
#define GRIDS_H_EUWJG

namespace zisa {

class TestGridFactory {
public:
  static std::string dbg() { return "grids/dbg.msh.h5"; }
  static std::string small() { return "grids/small.msh.h5"; }
  static std::string polytrope() { return "grids/unit_tests/polytrope.msh.h5"; }
  static std::string janka_polytrope() {
    return "grids/unit_tests/janka_polytrope.msh.h5";
  }

  static std::string unit_square(int i) {
    return string_format("grids/convergence/unit_square_%d/grid.msh.h5", i);
  }

  static std::string unit_square_with_halo(int i) {
    return string_format(
        "grids/convergence/unit_square_with_halo_%d/grid.msh.h5", i);
  }

  static std::vector<std::string> unit_square() {
    int n_grids = 4;
    auto pattern = [](int i) { return TestGridFactory::unit_square(i); };
    return generate_list(pattern, n_grids);
  }

  static std::string unit_cube(int i) {
    return string_format("grids/convergence/unit_cube_%d/grid.msh.h5", i);
  }

  static std::string unit_cube_with_halo(int i) {
    return string_format("grids/convergence/unit_cube_with_halo_%d/grid.msh.h5",
                         i);
  }

  static std::vector<std::string> unit_cube() {
    int n_grids = 3;
    auto pattern = [](int i) { return TestGridFactory::unit_cube(i); };
    return generate_list(pattern, n_grids);
  }

private:
  template <class F>
  static std::vector<std::string> generate_list(const F &f, int n_grids) {
    auto ret = std::vector<std::string>();
    ret.reserve(n_grids);

    for (int i = 0; i < n_grids; ++i) {
      ret.push_back(f(i));
    }

    return ret;
  }
};

} // namespace zisa

#endif /* end of include guard */
