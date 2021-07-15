// SPDX-License-Identifier: MIT
// Copyright (c) 2021 ETH Zurich, Luc Grosheintz-Laval

#include <zisa/parallelization/all_reduce.hpp>

namespace zisa {

double AllReduce::operator()(double local) const { return do_reduce(local); }

}