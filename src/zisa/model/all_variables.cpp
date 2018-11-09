/* Data-structure for all cell-centered variables exposed to Tyr.
 *
 * Authors: Luc Grosheintz <forbugrep@zoho.com>
 *    Date: 2016-06-13
 */
#include <zisa/model/all_variables.hpp>

namespace zisa {
AllVariables::AllVariables(const AllVariablesDimensions &dims) {
  allocate(dims);
}

void AllVariables::allocate(const AllVariablesDimensions &dims) {
  auto n_cells = dims.n_cells;

  cvars = GridVariables(shape_t<2>{n_cells, dims.n_cvars});
  avars = GridVariables(shape_t<2>{n_cells, dims.n_avars});
}

double AllVariables::operator[](int_t i) const {
  auto n_conserved_elements = cvars.size();

  if (i < n_conserved_elements) {
    return cvars[i];
  } else {
    return avars[i - n_conserved_elements];
  }
}

double &AllVariables::operator[](int_t i) {
  auto n_conserved_elements = cvars.size();

  if (i < n_conserved_elements) {
    return cvars[i];
  } else {
    return avars[i - n_conserved_elements];
  }
}

int_t AllVariables::size() const { return cvars.size() + avars.size(); }

AllVariablesDimensions AllVariables::dims(void) const {
  AllVariablesDimensions dims;

  dims.n_cells = cvars.shape(0);
  dims.n_cvars = cvars.shape(1);
  dims.n_avars = avars.shape(1);

  return dims;
}

void AllVariables::save(HDF5Writer &, const std::vector<std::string> &) const {
  LOG_ERR("Implement first.");
  // // conserved variables
  // cvars.split_save(writer, labels);

  // assert(avars.shape(1) == 0);
  // // advected variables
  // int n_avars = avars.shape(1);
  // writer.write_scalar(n_avars, "n_avars");

  // auto advected_labels = numbered_labels("mq%d", n_avars);
  // avars.split_save(writer, advected_labels);
}

// void AllVariables::load(HDF5Reader &reader,
//                         const std::vector<std::string> &labels) {
//   // conserved variables
//   cvars.split_load(reader, labels);

//   // advected variables
//   int n_avars =
//   reader.read_scalar<int>("n_avars");

//   auto advected_labels = numbered_labels("mq%d", n_avars);
//   avars.split_load(reader, advected_labels);
// }

} // namespace zisa
