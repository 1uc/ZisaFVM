#include <zisa/model/janka_eos.hpp>
namespace zisa {

std::string JankaEOSParams::str() const {
  return string_format("{%e, {%e, %e}, %e, %e}",
                       rho_bounce,
                       gamma[0],
                       gamma[1],
                       gamma_thermal,
                       E1);
}

void save(HDF5Writer &writer, const JankaEOSParams &params) {
  writer.write_scalar(params.rho_bounce, "rho_bounce");
  writer.write_scalar(params.gamma[0], "gamma1");
  writer.write_scalar(params.gamma[1], "gamma2");
  writer.write_scalar(params.gamma_thermal, "gamma_thermal");
  writer.write_scalar(params.E1, "E1");
}

JankaEOS make_janka_eos(const JankaEOSParams &params) {
  return JankaEOS(params);
}

JankaEOSParams make_default_janka_eos_params() {
  JankaEOSParams params;

  params.rho_bounce = 2e14;
  params.gamma = {1.33, 2.5};
  params.gamma_thermal = 1.5;
  params.E1 = 1.46925e15;

  return params;
}

JankaEOS make_default_janka_eos() {
  auto params = make_default_janka_eos_params();
  return make_janka_eos(params);
}

void save(HDF5Writer &writer, const JankaEOS &eos) {
  writer.write_string("JankaEOS", "eos");
  save(writer, eos.params());
}

std::string JankaEOS::str() const { return params().str(); }
}
