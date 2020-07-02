#ifndef EULER_FACTORY_H_BDRPI
#define EULER_FACTORY_H_BDRPI

#include <zisa/io/hdf5_serial_writer.hpp>
#include <zisa/model/euler_factory.hpp>

namespace zisa {
template <>
IdealGasEOS make_eos<IdealGasEOS>(const InputParameters &params) {
  auto eos_params = params["euler"]["eos"];

  double gamma = eos_params["gamma"];
  double R = eos_params["specific-gas-constant"];

  return IdealGasEOS(gamma, R);
}

JankaEOSParams make_janka_eos_params(const InputParameters &params) {
  const auto &eos_params = params["euler"]["eos"];
  LOG_ERR_IF(eos_params["mode"] != "janka", "Mismatching EOS.");

  JankaEOSParams janka_params;

  janka_params.rho_bounce = eos_params["rho_bounce"];
  janka_params.gamma[0] = eos_params["gamma1"];
  janka_params.gamma[1] = eos_params["gamma2"];
  janka_params.gamma_thermal = eos_params["gamma_thermal"];
  janka_params.E1 = eos_params["E1"];

  return janka_params;
}

template <>
JankaEOS make_eos<JankaEOS>(const InputParameters &params) {
  auto janka_params = make_janka_eos_params(params);
  return make_janka_eos(janka_params);
}

#if ZISA_HAS_HELMHOLTZ_EOS == 1
template <>
HelmholtzEOS make_eos<HelmholtzEOS>(const InputParameters &params) {
  const auto &eos_params = params["euler"]["eos"];
  auto table_path = std::string(eos_params["data"]);
  auto mass_mixing_ratio = std::vector<double>(eos_params["mass_mixing_ratio"]);
  auto mass_number = std::vector<double>(eos_params["mass_number"]);
  auto charge_number = std::vector<double>(eos_params["charge_number"]);

  return HelmholtzEOS(
      table_path, mass_mixing_ratio, mass_number, charge_number);
}
#endif

template <>
NoGravity make_gravity<NoGravity>(const InputParameters &params) {
  std::string mode = params["euler"]["gravity"]["mode"];
  LOG_ERR_IF(mode != "no_gravity", "Incompatible gravity.");

  return NoGravity{};
}

template <>
ConstantGravityRadial
make_gravity<ConstantGravityRadial>(const InputParameters &params) {
  LOG_ERR_IF(params["euler"]["gravity"]["mode"] != "constant",
             "Incompatible gravity.");

  double g = params["euler"]["gravity"]["g"];
  return ConstantGravityRadial(g);
}

template <>
PolytropeGravityRadial
make_gravity<PolytropeGravityRadial>(const InputParameters &input_params) {
  LOG_ERR_IF(input_params["euler"]["gravity"]["mode"] != "polytrope",
             "Incompatible gravity.");

  const auto &params = input_params["euler"]["gravity"];
  return PolytropeGravityRadial(params["rhoC"], params["K"], params["G"]);
}

template <>
PolytropeGravityWithJumpRadial make_gravity<PolytropeGravityWithJumpRadial>(
    const InputParameters &input_params) {

  LOG_ERR_IF(input_params["euler"]["gravity"]["mode"] != "polytrope_with_jump",
             "Incompatible gravity.");

  const auto &params = input_params["euler"]["gravity"];
  return PolytropeGravityWithJumpRadial(params["r_crit"],
                                        params["rhoC"],
                                        params["K_inner"],
                                        params["K_outer"],
                                        params["G"]);
}

template <>
RadialGravity make_gravity<RadialGravity>(const InputParameters &input_params) {
  LOG_ERR_IF(input_params["euler"]["gravity"]["mode"] != "radial_interpolation",
             "Incompatible gravity.");

  const auto &params = input_params["euler"]["gravity"];
  auto filename = params["profile"];

  auto reader = HDF5SerialReader(filename);
  auto points = array<double, 1>::load(reader, "x1");
  reader.open_group("auxiliary");
  auto grav_potential = array<double, 1>::load(reader, "grav._potential");

  return RadialGravity(std::move(points), std::move(grav_potential));
}

Euler<IdealGasEOS, ConstantGravityRadial> make_default_euler() {
  return Euler{IdealGasEOS{1.2, 2.3}, ConstantGravityRadial{0.99}};
}

} // namespace zisa

#endif /* end of include guard */
