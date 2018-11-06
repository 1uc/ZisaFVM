#ifndef CWENO_AO_H_CLLVU
#define CWENO_AO_H_CLLVU

#include <zisa/reconstruction/hybrid_weno.hpp>

namespace zisa {

class CWENO_AO : public HybridWENO {
private:
  using super = HybridWENO;

public:
  CWENO_AO(const std::shared_ptr<Grid> &grid,
           int_t i_cell,
           const HybridWENO_Params &params);

  auto reconstruct(const array<double, 1> &qbar) const -> decltype(hybridize());

  bool operator==(const CWENO_AO &other) const;
  bool operator!=(const CWENO_AO &other) const;

private:
  int_t k_high;
};

} // namespace zisa
#endif /* end of include guard */
