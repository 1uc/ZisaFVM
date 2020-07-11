#ifndef LOCAL_CFL_CONDITION_DECL_H_QY93E
#define LOCAL_CFL_CONDITION_DECL_H_QY93E

#include <zisa/grid/grid.hpp>
#include <zisa/model/all_variables.hpp>
#include <zisa/model/cfl_condition.hpp>

namespace zisa {

/// Compute the `dt_cfl` locally and then take the minimum.
template <class EOS>
class LocalCFL : public CFLCondition {
private:
  using super = CFLCondition;

protected:
  using cvars_t = euler_var_t;

public:
  LocalCFL(std::shared_ptr<Grid> grid,
           std::shared_ptr<Euler> model,
           std::shared_ptr<LocalEOSState<EOS>> local_eos,
           double cfl_number);

  virtual double operator()(const AllVariables &u) override;

protected:
  std::shared_ptr<Grid> grid;
  std::shared_ptr<Euler> euler;
  std::shared_ptr<LocalEOSState<EOS>> local_eos;
  double cfl_number;
};

} // namespace zisa

#endif
