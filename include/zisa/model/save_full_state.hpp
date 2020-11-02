#ifndef ZISA_SAVE_FULL_STATE_HPP
#define ZISA_SAVE_FULL_STATE_HPP

#include <zisa/io/hdf5.hpp>
#include <zisa/model/all_variables.hpp>
#include <zisa/model/euler.hpp>
#include <zisa/model/local_eos_state.hpp>

namespace zisa {

template <class EOS>
void save_full_state(HDF5Writer &writer,
                     const LocalEOSState<EOS> &local_eos,
                     const AllVariables &all_variables,
                     double t,
                     int_t n_steps) {

  ZISA_UNUSED(local_eos);

  auto labels = all_labels<typename Euler::cvars_t>();
  save_state(writer, all_variables, t, n_steps, labels);
  //  save_extended_state(writer, local_eos, all_variables);
}

template <class EOS>
void save_extended_state(HDF5Writer &writer,
                         const LocalEOSState<EOS> &local_eos,
                         const AllVariables &all_variables) {

  const auto &cvars = all_variables.cvars;
  const auto &avars = all_variables.avars;

  int_t n_cells = cvars.shape(0);
  int_t n_xvars = 7;

  // TODO this is not very efficient. It would be better to write to `n_xvars`
  //      separate arrays which could be serialized without unpacking.
  auto tmp_array = GridVariables(shape_t<2>{n_cells, n_xvars});

  for (int_t i = 0; i < n_cells; ++i) {
    auto u = euler_var_t(cvars(i));
    const auto &eos = *local_eos(i);
    auto full_xvars = eos.full_extra_variables(eos.rhoE(u));

    tmp_array(i, 0) = u(1) / u(0);
    tmp_array(i, 1) = u(2) / u(0);
    tmp_array(i, 2) = u(3) / u(0);
    tmp_array(i, 3) = full_xvars.p;
    tmp_array(i, 4) = full_xvars.a;
    tmp_array(i, 5) = full_xvars.h;
    tmp_array(i, 6) = full_xvars.s;
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

void save_extended_state(HDF5Writer &writer,
                         const JankaEOS &eos,
                         const AllVariables &all_variables);

}

#endif // ZISA_SAVE_FULL_STATE_HPP
