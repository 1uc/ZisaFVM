#include <zisa/parallelization/all_variables_scatterer.hpp>

namespace zisa {

AllVariablesScatterer::AllVariablesScatterer(
    std::unique_ptr<ArrayScatterer<double, 2>> cvars_scatterer,
    std::unique_ptr<ArrayScatterer<double, 2>> avars_scatterer,
    int_t n_cells_interior,
    int_t n_cells_with_halo)
    : cvars_scatterer(std::move(cvars_scatterer)),
      avars_scatterer(std::move(avars_scatterer)),
      n_cells_interior(n_cells_interior),
      n_cells_with_halo(n_cells_with_halo) {}

void AllVariablesScatterer::copy_local_patch(AllVariables &out,
                                             const AllVariables &in) {

  cvars_scatterer->copy_local_patch(subview(array_view(out.cvars)), in.cvars);

  if (in.avars.shape(1) != 0) {
    avars_scatterer->copy_local_patch(subview(array_view(out.avars)), in.avars);
  }
}

void AllVariablesScatterer::send(const AllVariables &all_variables) {
  cvars_scatterer->send(all_variables.cvars);

  if (all_variables.avars.shape(1) != 0) {
    avars_scatterer->send(all_variables.avars);
  }
}

void AllVariablesScatterer::receive(AllVariables &all_variables) {
  PRINT(all_variables.cvars.shape());
  cvars_scatterer->receive(subview(all_variables.cvars));

  if (all_variables.avars.shape(1) != 0) {
    avars_scatterer->receive(subview(all_variables.avars));
  }
}

bool AllVariablesScatterer::is_this_rank_scattering() const {
  bool is_scattering_cvars = cvars_scatterer->is_this_rank_scattering();

  if (avars_scatterer != nullptr) {
    bool is_scattering_avars = avars_scatterer->is_this_rank_scattering();
    LOG_ERR_IF(is_scattering_avars != is_scattering_cvars,
               "This needs to be implemented first.");
  }

  return is_scattering_cvars;
}

array_view<double, 2>
AllVariablesScatterer::subview(array_view<double, 2> view) const {
  auto shape = view.shape();
  shape[0] = n_cells_interior;

  return {shape, view.raw()};
}

AllVariables AllVariablesScatterer::scatter(const AllVariables &all_vars_full) {
  assert(is_this_rank_scattering());
  auto dims = all_vars_full.dims();

  auto all_vars_local
      = AllVariables({n_cells_with_halo, dims.n_cvars, dims.n_avars});

  send(all_vars_full);
  copy_local_patch(all_vars_local, all_vars_full);

  return all_vars_local;
}

AllVariables AllVariablesScatterer::scatter(int_t n_cvars, int_t n_avars) {
  assert(!is_this_rank_scattering());

  auto all_vars_local = AllVariables({n_cells_with_halo, n_cvars, n_avars});

  receive(all_vars_local);
  return all_vars_local;
}

}