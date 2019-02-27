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

static std::vector<std::string> numbered_labels(const std::string &pattern,
                                                int_t n_labels) {
  std::vector<std::string> labels;
  for (int_t i = 0; i < n_labels; ++i) {
    labels.push_back(string_format(pattern.c_str(), i));
  }

  return labels;
}

void save(HDF5Writer &writer,
          const AllVariables &all_variables,
          const std::vector<std::string> &labels) {

  const auto &cvars = all_variables.cvars;
  const auto &avars = all_variables.avars;

  // conserved variables
  save(writer, cvars, labels);

  // advected variables
  int_t n_avars = avars.shape(1);
  writer.write_scalar(n_avars, "n_avars");
  save(writer, avars, numbered_labels("mq%d", n_avars));
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
