#include <zisa/model/helmholtz_eos.hpp>

#if ZISA_HAS_HELMHOLTZ_EOS == 1

namespace zisa {
void initialize_helmholtz_eos(std::string &eos_table) {
  static auto done = [](std::string &eos_table) {
    int n_chars = int(eos_table.length());
    int status = -1;
    helmholtz_eos_readtable_c(eos_table.c_str(), &n_chars, &status);
    LOG_ERR_IF(
        status != 0,
        string_format("Reading EOS table failed. [%s]", eos_table.c_str()));

    return true;
  }(eos_table);

  ZISA_UNUSED(done);
}

void save(HDF5Writer & /* writer */, const HelmholtzEOS & /* eos */) {
  LOG_WARN("doing nothing.");
}
}

#endif
