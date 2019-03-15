#ifndef GLOBAL_RECONSTRUCTION_DECL_H_VHVV5
#define GLOBAL_RECONSTRUCTION_DECL_H_VHVV5

#include <zisa/config.hpp>
#include <zisa/grid/grid.hpp>
#include <zisa/memory/block_allocator.hpp>
#include <zisa/model/all_variables_fwd.hpp>
#include <zisa/model/euler_variables.hpp>
#include <zisa/reconstruction/hybrid_weno_params.hpp>
#include <zisa/reconstruction/local_reconstruction.hpp>
#include <zisa/reconstruction/weno_poly.hpp>

namespace zisa {

template <class CVars>
class GlobalReconstruction {
public:
  virtual ~GlobalReconstruction() = default;

  virtual void compute(const AllVariables &current_state) = 0;
  virtual CVars operator()(int_t i, const XYZ &x) const = 0;
};

template <class Equilibrium, class RC>
class EulerGlobalReconstruction : public GlobalReconstruction<euler_var_t> {
private:
  using cvars_t = euler_var_t;

public:
  EulerGlobalReconstruction(std::shared_ptr<Grid> grid,
                            const HybridWENOParams &params,
                            const Equilibrium &eq);

  const LocalReconstruction<Equilibrium, RC> &operator()(int_t i) const;
  virtual euler_var_t operator()(int_t i, const XYZ &x) const override;
  virtual void compute(const AllVariables &current_state) override;

  std::string str() const;

private:
  void set_qbar_local(array<cvars_t, 1> &qbar_local,
                      const AllVariables &current_state,
                      int_t i);

private:
  HybridWENOParams params;
  int_t max_stencil_size;

  array<LocalReconstruction<Equilibrium, RC>, 1> rc;
  std::shared_ptr<block_allocator<array<cvars_t, 1>>> allocator;
};

} // namespace zisa

#endif
