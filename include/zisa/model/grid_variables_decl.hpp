#ifndef GRID_VARIABLES_DECL_H_LKLNS
#define GRID_VARIABLES_DECL_H_LKLNS

#include <zisa/config.hpp>
#include <zisa/math/cartesian.hpp>
#include <zisa/memory/array.hpp>

namespace zisa {
class CVarExpr;
class CVarConstExpr;

class GridVariables : public array<double, 2> {
private:
  using super = array<double, 2>;

public:
  using super::super;

  inline CVarExpr operator()(int_t i);
  inline CVarConstExpr operator()(int_t i) const;

  inline double &operator()(int_t i, int_t k);
  inline double operator()(int_t i, int_t k) const;

  static GridVariables load(HDF5Reader &reader,
                            const std::vector<std::string> &labels);
};

void save(HDF5Writer &writer,
          const GridVariables &grid_variables,
          const std::vector<std::string> &labels);

} // namespace zisa

#endif /* end of include guard */
