// SPDX-License-Identifier: MIT
// Copyright (c) 2021 ETH Zurich, Luc Grosheintz-Laval

#ifndef CWENO_AO_H_CLLVU
#define CWENO_AO_H_CLLVU

#include <zisa/reconstruction/hybrid_weno.hpp>

namespace zisa {

class CWENO_AO : public HybridWENO {
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

private:
  template <class Poly>
  Poly reconstruct_impl(const array_view<double, 2, row_major> &rhs,
                        const array_view<Poly, 1> &polys,
                        const array_const_view<double, 2> &qbar) const;
};

} // namespace zisa
#endif /* end of include guard */
