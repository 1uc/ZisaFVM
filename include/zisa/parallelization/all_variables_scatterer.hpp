#ifndef ZISA_ALL_VARIABLES_SCATTERER_HPP_ZXZKQW
#define ZISA_ALL_VARIABLES_SCATTERER_HPP_ZXZKQW

#include <zisa/config.hpp>

#include <zisa/model/all_variables.hpp>
#include <zisa/parallelization/array_scatterer.hpp>

namespace zisa {

class AllVariablesScatterer {
public:
  AllVariablesScatterer(
      std::unique_ptr<ArrayScatterer<double, 2>> cvars_scatterer,
      std::unique_ptr<ArrayScatterer<double, 2>> avars_scatterer,
      int_t n_cells_local,
      int_t n_cells_global);

  AllVariables scatter(const AllVariables &all_vars_part);
  AllVariables scatter(int_t n_cvars, int_t n_avars);

  void copy_local_patch(AllVariables &out, const AllVariables &in);
  void send(const AllVariables &all_variables);
  void receive(AllVariables &all_variables);

  bool is_this_rank_scattering() const;

protected:
  array_view<double, 2, row_major>
  subview(array_view<double, 2, row_major> view) const;

private:
  std::unique_ptr<ArrayScatterer<double, 2>> cvars_scatterer;
  std::unique_ptr<ArrayScatterer<double, 2>> avars_scatterer;
  int_t n_cells_interior;
  int_t n_cells_with_halo;
};

}

#endif // ZISA_ALL_VARIABLES_SCATTERER_HPP
