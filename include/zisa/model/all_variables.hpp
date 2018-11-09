/* Array-like data-structures for FVM.
 *
 * Authors: Luc Grosheintz <forbugrep@zoho.com>
 *    Date: 2016-06-13
 */
#ifndef ALL_VARIABLES_H_NO06DDKF
#define ALL_VARIABLES_H_NO06DDKF

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

class AllVariables {
public:
  GridVariables cvars;  // conserved variables
  GridVariables avars;  // advected variables

public:
  AllVariables() = default;
  AllVariables(const AllVariablesDimensions &dims);

  double &operator[](int_t i);
  double operator[](int_t i) const;

  int_t size() const;

  AllVariablesDimensions dims(void) const;

  void save(HDF5Writer &writer, const std::vector<std::string> &labels) const;
  void load(HDF5Reader &reader, const std::vector<std::string> &labels);

protected:
  void allocate(const AllVariablesDimensions &dims);
};

} // namespace zisa

#endif /* end of include guard: ALL_VARIABLES_H_NO06DDKF */
