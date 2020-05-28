#ifndef ZISA_ALL_VARIABLES_GATHERER_HPP_IXKWO
#define ZISA_ALL_VARIABLES_GATHERER_HPP_IXKWO

#include <zisa/config.hpp>

#include <zisa/model/all_variables.hpp>
#include <zisa/parallelization/array_gatherer.hpp>

namespace zisa {

class AllVariablesGatherer {
public:
  AllVariablesGatherer(std::unique_ptr<ArrayGatherer<double, 2>> cvars_gatherer,
                       std::unique_ptr<ArrayGatherer<double, 2>> avars_gatherer,
                       int_t n_cells_local,
                       int_t n_cells_global);

  void copy_local_patch(AllVariables &out, const AllVariables &in);

  void send(const AllVariables &all_variables);

  void receive(AllVariables &all_variables);

  AllVariables gather(const AllVariables &all_vars_part);

  bool is_this_rank_gathering() const;

protected:
  array_const_view<double, 2, row_major>
  subview(const array_const_view<double, 2, row_major> &view) const;

private:
  std::unique_ptr<ArrayGatherer<double, 2>> cvars_gatherer;
  std::unique_ptr<ArrayGatherer<double, 2>> avars_gatherer;
  int_t n_cells_local;
  int_t n_cells_global;
};

}
#endif // ZISA_ALL_VARIABLES_GATHERER_HPP
