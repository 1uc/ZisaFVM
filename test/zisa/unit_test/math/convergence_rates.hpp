// SPDX-License-Identifier: MIT
// Copyright (c) 2021 ETH Zurich, Luc Grosheintz-Laval

/* Compute convergence rate.
 */

#ifndef BASIC_FUNCTIONS_H_ENWHK
#define BASIC_FUNCTIONS_H_ENWHK

#include <vector>

std::vector<double> convergence_rates(const std::vector<double> &resolution,
                                      const std::vector<double> &error);

#endif /* end of include guard */
