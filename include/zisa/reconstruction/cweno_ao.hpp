#ifndef CWENO_AO_H_CLLVU
#define CWENO_AO_H_CLLVU

#include <zisa/reconstruction/hybrid_weno.hpp>

namespace zisa {

class CWENO_AO : public HybridWENO {
private:
  using super = HybridWENO;

public:
  using super::super;

  auto reconstruct(const array<cvars_t, 1> &qbar) const -> decltype(hybridize());
};

} // namespace zisa
#endif /* end of include guard */
