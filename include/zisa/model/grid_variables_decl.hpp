// SPDX-License-Identifier: MIT
// Copyright (c) 2021 ETH Zurich, Luc Grosheintz-Laval

#ifndef GRID_VARIABLES_DECL_H_LKLNS
#define GRID_VARIABLES_DECL_H_LKLNS

#include "zisa/io/hierarchical_reader.hpp"
#include <zisa/config.hpp>
#include <zisa/loops/range.hpp>
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

  [[nodiscard]] static GridVariables
  load(HierarchicalReader &reader, const std::vector<std::string> &labels);

  static void load(HierarchicalReader &reader,
                   GridVariables &vars,
                   const std::vector<std::string> &labels);
};

void save(HierarchicalWriter &writer,
          const GridVariables &grid_variables,
          const std::vector<std::string> &labels);

class DereferenceConstGridVariables {
public:
  explicit DereferenceConstGridVariables(const GridVariables &grid_vars);

  inline CVarConstExpr item(int_t i) const;
  static constexpr bool has_item() { return true; }

private:
  const GridVariables &grid_vars;
};

Range<PlainIndexRange, DereferenceConstGridVariables>
cell_const_range(const GridVariables &grid_vars);

} // namespace zisa

#endif /* end of include guard */
