#ifndef LOCAL_CFL_CONDITION_DECL_H_QY93E
#define LOCAL_CFL_CONDITION_DECL_H_QY93E

#include <zisa/model/cfl_condition.hpp>
#include <zisa/model/all_variables.hpp>
#include <zisa/grid/grid.hpp>

namespace zisa {

/// Compute the `dt_cfl` locally and then take the minimum.
template <class Model>
class LocalCFL : public CFLCondition {
private:
  using super = CFLCondition;

protected:
  using cvars_t = typename Model::cvars_t;

public:
  LocalCFL(std::shared_ptr<Grid> grid, Model model, double cfl_number);

  virtual double operator()(const AllVariables &u) override;

protected:
  std::shared_ptr<Grid> grid;
  Model model;
  double cfl_number;
};

} // namespace zisa

#endif
