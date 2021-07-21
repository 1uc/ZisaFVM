// SPDX-License-Identifier: MIT
// Copyright (c) 2021 ETH Zurich, Luc Grosheintz-Laval

//
// Created by lucg on 11/17/20.
//

#ifndef ZISA_LOAD_FULL_STATE_HPP_BWCBU
#define ZISA_LOAD_FULL_STATE_HPP_BWCBU

#include <zisa/model/models.hpp>
#include <zisa/io/hdf5_writer.hpp>
#include <zisa/model/all_variables.hpp>
#include <zisa/model/euler.hpp>

namespace zisa {

std::pair<double, int_t> load_full_state(HDF5Reader &reader,
                     AllVariables &all_variables);

}
#endif // ZISA_LOAD_FULL_STATE_HPP
