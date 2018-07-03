/* Data-structure for all cell-centered variables exposed to Tyr.
 *
 * Authors: Luc Grosheintz <forbugrep@zoho.com>
 *    Date: 2016-06-13
 */
#ifndef ALL_VARIABLES_CPP_WLARDDKM
#define ALL_VARIABLES_CPP_WLARDDKM

#include "tyr/model/all_variables.h"
#include "esp/config.h"
#include "esp/grid/structured_grid.h"
#include "esp/memory/fill.h"

namespace tyr {
void GridVariables::split_save(HDF5Writer &writer,
                               const std::vector<std::string> &labels) const {
  Array<double, 1> component{shape[0] * shape[1]};
  for (int k = 0; k < int(labels.size()); ++k) {
    for (int l = 0; l < shape[0]; ++l) {
      for (int i = 0; i < shape[1]; ++i) {
        component[linear_index(l, i)] = (*this)(l, i, k);
      }
    }

    component.save(writer, labels.at(k));
  }
  component.free();
}

void GridVariables::split_save(HDF5Writer &writer,
                               const StructuredGrid &grid,
                               const std::vector<std::string> &labels) const {
  Array<double, 3> component(grid.shape);
  for (int k = 0; k < int(labels.size()); ++k) {
    for (auto z_range = grid.range<2>(0); !z_range.is_finished(); ++z_range) {
      for (auto y_range = grid.range<1>(0); !y_range.is_finished(); ++y_range) {
        for (auto x_range = grid.range<0>(0); !x_range.is_finished();
             ++x_range) {
          int ix = x_range.idx, iy = y_range.idx, iz = z_range.idx;
          component(ix, iy, iz) = (*this)(grid.cell_index(ix, iy, iz), k);
        }
      }
    }

    component.save(writer, labels.at(k));
  }
  component.free();
}

void GridVariables::split_load(HDF5Reader &reader,
                               const std::vector<std::string> &labels) {
  Array<double, 1> component;

  for (int k = 0; k < int(labels.size()); ++k) {
    component.load(reader, labels.at(k));

    for (int l = 0; l < this->shape[0]; ++l) {
      for (int i = 0; i < this->shape[1]; ++i) {
        (*this)(l, i, k) = component[linear_index(l, i)];
      }
    }
  }

  component.free();
}

void GridVariables::fill_host(double value) {
  tyr::fill(CStyleMemoryManager(), this->begin(), this->end(), value);
}

std::ostream &operator<<(std::ostream &os, const AllVariablesDimensions &dims) {
  return os << dims.str();
}

AllVariables::AllVariables(const MemoryManager &mm,
                           const AllVariablesDimensions &dims) {
  allocate(mm, dims);
}

double AllVariables::operator[](int i) const {
  int n_conserved_elements = conserved_variables.size();

  if (i < n_conserved_elements) {
    return conserved_variables[i];
  } else {
    return advected_variables[i - n_conserved_elements];
  }
}

double &AllVariables::operator[](int i) {
  int n_conserved_elements = conserved_variables.size();

  if (i < n_conserved_elements) {
    return conserved_variables[i];
  } else {
    return advected_variables[i - n_conserved_elements];
  }
}

int AllVariables::size() const {
  return conserved_variables.size() + advected_variables.size();
}

void AllVariables::allocate(const MemoryManager &mm,
                            const AllVariablesDimensions &dims) {
  int n_layers = dims.n_layers;
  int n_cells = dims.n_cells;

  conserved_variables.allocate(mm,
                               {n_layers, n_cells, dims.n_conserved_variables});
  advected_variables.allocate(mm,
                              {n_layers, n_cells, dims.n_advected_variables});
}

void AllVariables::free(const MemoryManager &mm) {
  conserved_variables.free(mm);
  advected_variables.free(mm);
}

AllVariablesDimensions AllVariables::dims(void) const {
  AllVariablesDimensions dims;

  dims.n_layers = conserved_variables.shape[0];
  dims.n_cells = conserved_variables.shape[1];

  dims.n_conserved_variables = conserved_variables.shape[2];
  dims.n_advected_variables = advected_variables.shape[2];

  return dims;
}

void AllVariables::fill_host(double value) {
  conserved_variables.fill_host(value);
  advected_variables.fill_host(value);
}

void AllVariables::save(HDF5Writer &writer,
                        const std::vector<std::string> &labels) const {
  // conserved variables
  conserved_variables.split_save(writer, labels);

  // advected variables
  int n_advected_variables = advected_variables.shape[2];
  writer.write_scalar(n_advected_variables, "n_advected_variables");

  auto advected_labels = numbered_labels("mq%d", n_advected_variables);
  advected_variables.split_save(writer, advected_labels);
}

void AllVariables::save(HDF5Writer &writer,
                        const StructuredGrid &grid,
                        const std::vector<std::string> &labels) const {
  // conserved variables
  conserved_variables.split_save(writer, grid, labels);

  // advected variables
  int n_advected_variables = advected_variables.shape[2];
  writer.write_scalar(n_advected_variables, "n_advected_variables");

  auto advected_labels = numbered_labels("mq%d", n_advected_variables);
  advected_variables.split_save(writer, grid, advected_labels);
}

void AllVariables::load(HDF5Reader &reader,
                        const std::vector<std::string> &labels) {
  // conserved variables
  conserved_variables.split_load(reader, labels);

  // advected variables
  int n_advected_variables = reader.read_scalar<int>("n_advected_variables");

  auto advected_labels = numbered_labels("mq%d", n_advected_variables);
  advected_variables.split_load(reader, advected_labels);
}

void copy_from_host(AllVariables &dest,
                    const AllVariables &src,
                    Accelerator accelerator) {
  copy_from_host(
      dest.conserved_variables, src.conserved_variables, accelerator);
  copy_from_host(dest.advected_variables, src.advected_variables, accelerator);
}

void copy_to_host(AllVariables &dest, const AllVariables &src) {
  copy_to_host(dest.conserved_variables, src.conserved_variables);
  copy_to_host(dest.advected_variables, src.advected_variables);
}

void copy_to_device(AllVariables &dest, const AllVariables &src) {
  copy_to_device(dest.conserved_variables, src.conserved_variables);
  copy_to_device(dest.advected_variables, src.advected_variables);
}

void copy_device_to_device(AllVariables &dest, const AllVariables &src) {
  copy_device_to_device(dest.conserved_variables, src.conserved_variables);
  copy_device_to_device(dest.advected_variables, src.advected_variables);
}

} // namespace tyr
#endif /* end of include guard: ALL_VARIABLES_CPP_WLARDDKM */
