/* List of grids that can be used for convergence studies.
 */

#ifndef GRIDS_H_EUWJG
#define GRIDS_H_EUWJG

namespace zisa {

class TestGridFactory {

public:
  static std::vector<std::string> unit_square() {
    int n_grids = 4;

    auto ret = std::vector<std::string>();
    ret.reserve(n_grids);

    for (int i = 0; i < n_grids; ++i) {
      ret.push_back(string_format("grids/convergence/unit_square_%d.msh", i));
    }

    return ret;
  }
};

} // namespace zisa

#endif /* end of include guard */
