// SPDX-License-Identifier: MIT
// Copyright (c) 2021 ETH Zurich, Luc Grosheintz-Laval

#include <zisa/model/ideal_gas_eos.hpp>

namespace zisa {

void save(HierarchicalWriter &writer, const IdealGasEOS &eos) {

  writer.open_group("eos");
  writer.write_scalar(eos.gamma(), "gamma");
  writer.write_scalar(eos.specific_gas_constant(), "specific_gas_constant");
  writer.close_group();

}

} // zisa
