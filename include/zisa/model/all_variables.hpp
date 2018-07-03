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

namespace zisa {
class GridVariables : public array<double, 2> {
private:
  using super = array<double, 2>;

public:
  using super::super;

  using super::operator();

  template <class Variable>
  ANY_DEVICE_INLINE Variable read(Variable &var, int i) const {
    for (int k = 0; k < (*this).shape(1); k++) {
      var[k] = (*this)(i, k);
    }
  }

  template <class Variable>
  ANY_DEVICE_INLINE Variable read(int i) const {
    Variable var;
    read(var, i);
    return var;
  }

  template <class Variable>
  ANY_DEVICE_INLINE void write(const Variable &var, int i) {
    for (int k = 0; k < (*this).shape(1); k++) {
      (*this)(i, k) = var[k];
    }
  }

  /// Save each physical variable separately.
  /** The array is split along the last dimension and stored as separate
   *  arrays in the HDF5 file. The name of each array is given in
   *  `labels`.
   *
   *  @note This stores the array as an 1D structure. It's suitable for
   *    unstructured grids.
   */
  void split_save(HDF5Writer &writer,
                  const std::vector<std::string> &labels) const;

  /// Load components stored in separate HDF5 arrays.
  /** Load one array for every component and fuse them into one array.
   */
  void split_load(HDF5Reader &reader, const std::vector<std::string> &labels);

  /// Fill the object with `values`.
  void fill_host(double value);
};

struct AllVariablesDimensions {
  int n_cells;
  int n_conserved_variables;
  int n_advected_variables;
};

std::ostream &operator<<(std::ostream &os, const AllVariablesDimensions &dims);

class AllVariables {
public:
  GridVariables conserved_variables;
  GridVariables advected_variables;

public:
  AllVariables() = default;
  AllVariables(const device_type &device, const AllVariablesDimensions &dims);

  double &operator[](int i);
  double operator[](int i) const;

  void allocate(const device_type &device, const AllVariablesDimensions &dims);

  int_t size() const;

  AllVariablesDimensions dims(void) const;

  /// Fill the object with `values`.
  void fill(double value);

  void save(HDF5Writer &writer, const std::vector<std::string> &labels) const;
  void load(HDF5Reader &reader, const std::vector<std::string> &labels);
};

} // namespace zisa

#endif /* end of include guard: ALL_VARIABLES_H_NO06DDKF */
