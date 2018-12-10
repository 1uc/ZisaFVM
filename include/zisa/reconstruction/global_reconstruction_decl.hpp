#ifndef GLOBAL_RECONSTRUCTION_DECL_H_VHVV5
#define GLOBAL_RECONSTRUCTION_DECL_H_VHVV5

#include <zisa/config.hpp>

#include <zisa/grid/grid.hpp>
#include <zisa/reconstruction/hybrid_weno_params.hpp>
#include <zisa/reconstruction/weno_poly.hpp>
#include <zisa/model/all_variables_fwd.hpp>

namespace zisa {

template <class RC>
class GlobalReconstruction {
public:
  GlobalReconstruction(std::shared_ptr<Grid> grid,
                       const HybridWENO_Params &params,
                       int_t n_vars);

  const WENOPoly &operator()(int_t i, int_t k) const;
  void compute(const AllVariables &current_state);

  std::string str() const;

private:
  void set_qbar_local(const AllVariables &current_state, int_t i, int_t k);

private:
  HybridWENO_Params params;

  array<RC, 1> rc;
  array<WENOPoly, 2> polys;

  mutable array<double, 1> qbar_local;
};

} // namespace zisa

#endif
