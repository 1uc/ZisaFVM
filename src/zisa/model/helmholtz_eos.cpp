#include <zisa/model/helmholtz_eos.hpp>

#if ZISA_HAS_HELMHOLTZ_EOS == 1

void zisa::save(HDF5Writer & /* writer */, const HelmholtzEOS & /* eos */) {
  LOG_WARN("doing nothing.");
}

#endif
