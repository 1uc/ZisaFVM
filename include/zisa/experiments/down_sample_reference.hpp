// SPDX-License-Identifier: MIT
// Copyright (c) 2021 ETH Zurich, Luc Grosheintz-Laval

#ifndef ZISA_DOWN_SAMPLE_REFERENCE_HPP
#define ZISA_DOWN_SAMPLE_REFERENCE_HPP
#include <zisa/config.hpp>

#include <zisa/cli/input_parameters.hpp>
#include <zisa/math/reference_solution.hpp>

namespace zisa {

void down_sample_euler_reference(
    const ReferenceSolution &reference_solution,
    const std::vector<std::string> &coarse_grids,
    const std::function<std::shared_ptr<Grid>(const std::string &, int_t)>
        &factory,
    const std::string &filename);

}

#endif // ZISA_DOWN_SAMPLE_REFERENCE_HPP
