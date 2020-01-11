#include <zisa/parallelization/all_variables_gatherer.hpp>

namespace zisa {

AllVariablesGatherer::AllVariablesGatherer(
    std::unique_ptr<ArrayGatherer<double, 2>> cvars_gatherer,
    std::unique_ptr<ArrayGatherer<double, 2>> avars_gatherer,
    int_t n_cells_local)
    : cvars_gatherer(std::move(cvars_gatherer)),
      avars_gatherer(std::move(avars_gatherer)),
      n_cells_local(n_cells_local) {}

void AllVariablesGatherer::copy_local_patch(AllVariables &out,
                                            const AllVariables &in) {

  cvars_gatherer->copy_local_patch(out.cvars, subview(in.cvars));

  if (in.avars.shape(1) != 0) {
    avars_gatherer->copy_local_patch(out.avars, subview(in.avars));
  }
}

void AllVariablesGatherer::send(const AllVariables &all_variables) {
  cvars_gatherer->send(subview(all_variables.cvars));

  if (all_variables.avars.shape(1) != 0) {
    avars_gatherer->send(subview(all_variables.avars));
  }
}

void AllVariablesGatherer::receive(AllVariables &all_variables) {
  cvars_gatherer->receive(all_variables.cvars);

  if (all_variables.avars.shape(1) != 0) {
    avars_gatherer->receive(all_variables.avars);
  }
}

bool AllVariablesGatherer::is_this_rank_gathering() const {
  bool is_gathering_cvars = cvars_gatherer->is_this_rank_gathering();

  if (avars_gatherer != nullptr) {
    bool is_gathering_avars = avars_gatherer->is_this_rank_gathering();
    LOG_ERR_IF(is_gathering_avars != is_gathering_cvars,
               "This needs to be implemented first.");
  }

  return is_gathering_cvars;
}

array_const_view<double, 2>
AllVariablesGatherer::subview(const array_const_view<double, 2> &view) const {
  auto shape = view.shape();
  shape[0] = n_cells_local;

  return {shape, view.raw()};
}
}