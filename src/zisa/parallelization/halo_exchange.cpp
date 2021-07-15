// SPDX-License-Identifier: MIT
// Copyright (c) 2021 ETH Zurich, Luc Grosheintz-Laval

#include <zisa/parallelization/halo_exchange.hpp>

namespace zisa {
void NoHaloExchange::operator()(AllVariables &) { return; }
void NoHaloExchange::wait() { return; }
}
