#ifndef GLOBAL_RECONSTRUCTION_DECL_H_VHVV5
#define GLOBAL_RECONSTRUCTION_DECL_H_VHVV5

#include <zisa/config.hpp>

#include <zisa/grid/grid.hpp>
#include <zisa/model/all_variables_fwd.hpp>
#include <zisa/model/euler_variables.hpp>
#include <zisa/reconstruction/hybrid_weno_params.hpp>
#include <zisa/reconstruction/local_reconstruction.hpp>
#include <zisa/reconstruction/weno_poly.hpp>

namespace zisa {

template <class Equilibrium, class RC>
class GlobalReconstruction {
private:
  using cvars_t = euler_var_t;

public:
  GlobalReconstruction(std::shared_ptr<Grid> grid,
                       const HybridWENO_Params &params,
                       const Equilibrium &eq);

  const LocalReconstruction<Equilibrium, RC> &operator()(int_t i) const;
  void compute(const AllVariables &current_state);

  std::string str() const;

private:
  void set_qbar_local(const AllVariables &current_state, int_t i);

private:
  HybridWENO_Params params;

  array<LocalReconstruction<Equilibrium, RC>, 1> rc;
  mutable std::vector<array<cvars_t, 1>> qbar_local;
};

} // namespace zisa

#endif
