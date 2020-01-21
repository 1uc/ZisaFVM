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
  using super::super;

  auto reconstruct(array<double, 2, row_major> &rhs,
                   array<WENOPoly, 1> &polys,
                   const array<cvars_t, 1> &qbar) const
      -> decltype(hybridize(polys));
};

} // namespace zisa
#endif /* end of include guard */
