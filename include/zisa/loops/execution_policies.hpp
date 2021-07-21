// SPDX-License-Identifier: MIT
// Copyright (c) 2021 ETH Zurich, Luc Grosheintz-Laval

#ifndef ZISA_EXECUTION_POLICIES_HPP
#define ZISA_EXECUTION_POLICIES_HPP

namespace zisa {

struct omp_policy {};
struct serial_policy {};

#if ZISA_HAS_OPENMP == 1
using default_execution_policy = omp_policy;
#else
using default_execution_policy = serial_policy;
#endif

}

#endif // ZISA_EXECUTION_POLICIES_HPP
