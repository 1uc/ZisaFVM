// SPDX-License-Identifier: MIT
// Copyright (c) 2021 ETH Zurich, Luc Grosheintz-Laval

#ifndef ZISA_LOCAL_EOS_STATE_HPP
#define ZISA_LOCAL_EOS_STATE_HPP

#include <zisa/model/all_variables.hpp>
#include <zisa/model/equation_of_state.hpp>

namespace zisa {

template <class EOS>
class LocalEOSState;

template <>
class LocalEOSState<IdealGasEOS> {
public:
  using eos_t = IdealGasEOS;

public:
  LocalEOSState(double gamma, double specific_gas_constant)
      : eos(std::make_shared<IdealGasEOS>(gamma, specific_gas_constant)) {}

  void compute(const AllVariables &) {
    // nothing to do.
  }

  void compute(const AllVariables &, const array_const_view<int_t, 1> &) {
    // nothing to do.
  }

  std::shared_ptr<IdealGasEOS> operator()(int_t /* i */) const { return eos; }

private:
  std::shared_ptr<IdealGasEOS> eos;
};

#if ZISA_HAS_HELMHOLTZ_EOS == 1
template <>
class LocalEOSState<HelmholtzEOS> {
public:
  using eos_t = HelmholtzEOS;

public:
  LocalEOSState() = delete;
  LocalEOSState(int_t n_cells,
                std::vector<double> mass_number,
                std::vector<double> charge_number)
      : local_eos(shape_t<1>(n_cells)),
        mass_fraction(mass_number.size()),
        mass_number(std::move(mass_number)),
        charge_number(std::move(charge_number)) {

    for (int_t i = 0; i < n_cells; ++i) {
      local_eos[i] = std::make_shared<HelmholtzEOS>();
    }
  }

  void compute(const AllVariables &all_variables) {
    const auto &avars = all_variables.avars;
    const auto &cvars = all_variables.cvars;

    auto n_cells = avars.shape(0);
    auto n_avars = avars.shape(1);

    LOG_ERR_IF(n_cells != local_eos.shape(0),
               string_format(
                   "Sizes don't match: %d != %d", n_cells, local_eos.shape(0)));

    for (int_t i = 0; i < n_cells; ++i) {
      for (int_t k = 0; k < n_avars; ++k) {
        mass_fraction[k] = zisa::max(0.0, avars(i, k) / cvars(i, 0));
      }

      *local_eos[i] = HelmholtzEOS(mass_fraction, mass_number, charge_number);
    }
  }

  void compute(const AllVariables &all_variables,
               const array_const_view<int_t, 1> &cells) {

    const auto &avars = all_variables.avars;
    const auto &cvars = all_variables.cvars;

    auto n_cells = cells.size();
    auto n_avars = avars.shape(1);

    LOG_ERR_IF(avars.shape(0) != local_eos.shape(0),
               string_format(
                   "Sizes don't match: %d != %d", n_cells, local_eos.shape(0)));

    // TODO not OpenMP parallel
    for (int_t ii = 0; ii < n_cells; ++ii) {
      auto i = cells[ii];
      for (int_t k = 0; k < n_avars; ++k) {
        mass_fraction[k] = zisa::max(0.0, avars(i, k) / cvars(i, 0));
      }

      *local_eos[i] = HelmholtzEOS(mass_fraction, mass_number, charge_number);
    }
  }

  std::shared_ptr<HelmholtzEOS> operator()(int_t i) const {
    return local_eos[i];
  }

private:
  array<std::shared_ptr<HelmholtzEOS>, 1> local_eos;
  std::vector<double> mass_fraction;
  std::vector<double> mass_number;
  std::vector<double> charge_number;
};
#endif

}

#endif // ZISA_LOCAL_EOS_STATE_HPP
