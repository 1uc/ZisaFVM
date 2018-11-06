/*
 *
 */

#ifndef WENO_AO_H_A1UM2
#define WENO_AO_H_A1UM2

#include <zisa/config.hpp>

#include <zisa/reconstruction/hybrid_weno.hpp>
#include <zisa/reconstruction/hybrid_weno_params.hpp>

namespace zisa {

class WENO_AO : public HybridWENO {
private:
  using super = HybridWENO;

public:
  WENO_AO(const std::shared_ptr<Grid> &grid,
          int_t i_cell,
          const HybridWENO_Params &params);

  auto reconstruct(const array<double, 1> &qbar) const -> decltype(hybridize());
};

} // namespace zisa
#endif /* end of include guard */
