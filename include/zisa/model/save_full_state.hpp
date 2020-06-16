#ifndef ZISA_SAVE_FULL_STATE_HPP
#define ZISA_SAVE_FULL_STATE_HPP

#include <zisa/io/hdf5.hpp>
#include <zisa/model/all_variables.hpp>
#include <zisa/model/euler.hpp>

namespace zisa {

template <class EOS, class Gravity>
void save_full_state(HDF5Writer &writer,
                     const Euler<EOS, Gravity> &euler,
                     const AllVariables &all_variables,
                     double t,
                     int_t n_steps) {

  save_state(writer, euler, all_variables, t, n_steps);
  save_extended_state(writer, euler, all_variables);
}

template <class EOS, class Gravity>
void save_extended_state(HDF5Writer &writer,
                         const Euler<EOS, Gravity> &euler,
                         const AllVariables &all_variables) {

  const auto &cvars = all_variables.cvars;
  const auto &avars = all_variables.avars;
  const auto &eos = euler.eos;

  int_t n_cells = cvars.shape(0);
  int_t n_xvars = 7;

  // TODO this is not very efficient. It would be better to write to `n_xvars`
  //      separate arrays which could be serialized without unpacking.
  auto tmp_array = GridVariables(shape_t<2>{n_cells, n_xvars});

  for (int_t i = 0; i < n_cells; ++i) {
    auto u = euler_var_t(cvars(i));
    auto rhoP = eos.rhoP(u);
    auto cs = eos.sound_speed(rhoP);
    auto theta = eos.enthalpy_entropy(rhoP);

    tmp_array(i, 0) = u(1) / u(0);
    tmp_array(i, 1) = u(2) / u(0);
    tmp_array(i, 2) = u(3) / u(0);
    tmp_array(i, 3) = rhoP.p();
    tmp_array(i, 4) = cs;
    tmp_array(i, 5) = theta.h();
    tmp_array(i, 6) = theta.s();
  }

  auto labels = std::vector<std::string>{"v1", "v2", "v3", "p", "cs", "h", "s"};
  save(writer, tmp_array, labels);

  auto avars_array = GridVariables(avars.shape());
  auto n_avars = avars_array.shape(1);

  for (int_t i = 0; i < n_cells; ++i) {
    for (int_t k_var = 0; k_var < n_avars; ++k_var) {
      avars_array(i, k_var) = avars(i, k_var) / cvars(i, 0);
    }
  }

  auto avar_labels = numbered_labels("q%d", n_avars);
  save(writer, avars_array, avar_labels);
}

template <class Gravity>
void save_extended_state(HDF5Writer &writer,
                         const Euler<JankaEOS, Gravity> &euler,
                         const AllVariables &all_variables) {

  const auto &cvars = all_variables.cvars;
  const auto &eos = euler.eos;

  int_t n_cells = cvars.shape(0);
  int_t n_xvars = 8;

  auto tmp_array = GridVariables(shape_t<2>{n_cells, n_xvars});

  for (int_t i = 0; i < n_cells; ++i) {
    auto u = euler_var_t(cvars(i));
    auto rhoP = eos.rhoP(u);
    auto cs = eos.sound_speed(rhoP);
    auto theta = eos.enthalpy_entropy(rhoP);

    tmp_array(i, 0) = rhoP.p();
    tmp_array(i, 1) = cs;
    tmp_array(i, 2) = theta.h();
    tmp_array(i, 3) = theta.s();
    tmp_array(i, 4) = eos.polytropic_energy(u[0]);
    tmp_array(i, 5) = eos.thermal_energy(rhoP);
    tmp_array(i, 6) = eos.polytropic_pressure(u[0]);
    tmp_array(i, 7) = eos.thermal_pressure(rhoP);
  }

  auto labels = std::vector<std::string>{
      "p", "cs", "h", "s", "E_p", "E_th", "p_p", "p_th"};
  save(writer, tmp_array, labels);
}
}

#endif // ZISA_SAVE_FULL_STATE_HPP
