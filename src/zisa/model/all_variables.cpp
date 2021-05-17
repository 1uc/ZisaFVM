#include <zisa/model/all_variables.hpp>

namespace zisa {
bool operator==(const AllVariablesDimensions &a,
                const AllVariablesDimensions &b) {
  return (a.n_cells == b.n_cells) && (a.n_cvars == b.n_cvars)
         && (a.n_avars == b.n_avars);
}

std::ostream &operator<<(std::ostream &os, const AllVariablesDimensions &dims) {
  os << "{ " << dims.n_cells << ", " << dims.n_cvars << ", " << dims.n_avars
     << "}";
  return os;
}
AllVariables::AllVariables(GridVariables cvars, GridVariables avars)
    : cvars(std::move(cvars)), avars(std::move(avars)) {}

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

AllVariablesDimensions AllVariables::dims() const {
  AllVariablesDimensions dims{};

  dims.n_cells = cvars.shape(0);
  dims.n_cvars = cvars.shape(1);
  dims.n_avars = avars.shape(1);

  return dims;
}

std::vector<std::string> numbered_labels(const std::string &pattern,
                                         int_t n_labels) {
  std::vector<std::string> labels;
  for (int_t i = 0; i < n_labels; ++i) {
    labels.push_back(string_format(pattern.c_str(), i));
  }

  return labels;
}

void save(HierarchicalWriter &writer,
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

[[nodiscard]] AllVariables
AllVariables::load(HierarchicalReader &reader, const std::vector<std::string> &labels) {

  auto all_vars = AllVariables{};
  all_vars.cvars = GridVariables::load(reader, labels);

  auto n_avars = reader.read_scalar<int_t>("n_avars");
  auto avar_labels = numbered_labels("mq%d", n_avars);
  all_vars.avars = GridVariables::load(reader, avar_labels);

  return all_vars;
}

void AllVariables::load(HierarchicalReader &reader,
                        AllVariables &all_vars,
                        const std::vector<std::string> &labels) {

  GridVariables::load(reader, all_vars.cvars, labels);

  auto n_avars = reader.read_scalar<int_t>("n_avars");
  auto avar_labels = numbered_labels("mq%d", n_avars);
  GridVariables::load(reader, all_vars.avars, avar_labels);
}

} // namespace zisa
