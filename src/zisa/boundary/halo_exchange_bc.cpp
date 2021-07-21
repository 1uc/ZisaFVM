// SPDX-License-Identifier: MIT
// Copyright (c) 2021 ETH Zurich, Luc Grosheintz-Laval

#include <zisa/boundary/halo_exchange_bc.hpp>

namespace zisa {

HaloExchangeBC::HaloExchangeBC(std::shared_ptr<BoundaryCondition> local_bc,
                               std::shared_ptr<HaloExchange> halo_exchange)
    : local_bc(std::move(local_bc)), halo_exchange(std::move(halo_exchange)) {}

void HaloExchangeBC::apply(AllVariables &u, double t) {
  local_bc->apply(u, t);
  (*halo_exchange)(u);
}

std::string HaloExchangeBC::str() const {
  return "HaloExchangeBC with:\n" + indent_block(1, local_bc->str());
}

}