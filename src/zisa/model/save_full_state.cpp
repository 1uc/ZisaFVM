#include <zisa/model/save_full_state.hpp>

namespace zisa {


void save_extended_state(HierarchicalWriter &writer,
                         const JankaEOS &eos,
                         const AllVariables &all_variables) {

  const auto &cvars = all_variables.cvars;

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
