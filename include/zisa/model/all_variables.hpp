/* Array-like data-structures for FVM.
 *
 * Authors: Luc Grosheintz <forbugrep@zoho.com>
 *    Date: 2016-06-13
 */
#ifndef ALL_VARIABLES_H_NO06DDKF
#define ALL_VARIABLES_H_NO06DDKF

#include "zisa/io/hierarchical_reader.hpp"
#include <zisa/config.hpp>
#include <zisa/math/cartesian.hpp>
#include <zisa/memory/array.hpp>
#include <zisa/model/grid_variables.hpp>

namespace zisa {

struct AllVariablesDimensions {
  int_t n_cells;
  int_t n_cvars;
  int_t n_avars;
};

std::ostream &operator<<(std::ostream &os, const AllVariablesDimensions &dims);

bool operator==(const AllVariablesDimensions &a,
                const AllVariablesDimensions &b);

inline bool operator!=(const AllVariablesDimensions &a,
                       const AllVariablesDimensions &b) {
  return !(a == b);
}

class AllVariables {
public:
  GridVariables cvars; // conserved variables
  GridVariables avars; // advected variables

public:
  AllVariables() = default;
  AllVariables(GridVariables cvars, GridVariables avars);
  AllVariables(const AllVariablesDimensions &dims);

  double &operator[](int_t i);
  double operator[](int_t i) const;

  int_t size() const;
  AllVariablesDimensions dims() const;

  [[nodiscard]] static AllVariables
  load(HierarchicalReader &reader, const std::vector<std::string> &labels);

  static void load(HierarchicalReader &reader,
                   AllVariables &all_vars,
                   const std::vector<std::string> &labels);

protected:
  void allocate(const AllVariablesDimensions &dims);
};

void save(HierarchicalWriter &writer,
          const AllVariables &all_variables,
          const std::vector<std::string> &labels);

std::vector<std::string> numbered_labels(const std::string &pattern,
                                         int_t n_labels);

} // namespace zisa

#endif /* end of include guard: ALL_VARIABLES_H_NO06DDKF */
