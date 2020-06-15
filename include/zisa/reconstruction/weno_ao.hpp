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

  WENOPoly reconstruct(const array_view<double, 2, row_major> &rhs,
                       const array_view<WENOPoly, 1> &polys,
                       const array_const_view<cvars_t, 1> &qbar) const;

  ScalarPoly reconstruct(const array_view<double, 2, row_major> &rhs,
                         const array_view<ScalarPoly, 1> &polys,
                         const array_const_view<double, 1> &qbar) const;
};

} // namespace zisa
#endif /* end of include guard */
