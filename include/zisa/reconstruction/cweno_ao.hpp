#ifndef CWENO_AO_H_CLLVU
#define CWENO_AO_H_CLLVU

#include <zisa/reconstruction/hybrid_weno.hpp>

namespace zisa {

class CWENO_AO : public HybridWENO {
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
