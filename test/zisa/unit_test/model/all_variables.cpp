// SPDX-License-Identifier: MIT
// Copyright (c) 2021 ETH Zurich, Luc Grosheintz-Laval

#include <random>

#include <zisa/io/hdf5_serial_writer.hpp>
#include <zisa/model/all_variables.hpp>
#include <zisa/model/euler.hpp>
#include <zisa/model/euler_factory.hpp>
#include <zisa/testing/testing_framework.hpp>

TEST_CASE("GridVariables", "[memory]") {
  zisa::int_t n_cells = 10;
  zisa::int_t n_vars = 5;

  auto u = zisa::GridVariables(zisa::shape_t<2>{n_cells, n_vars});
  std::fill(u.begin(), u.end(), 0.0);

  auto ua = zisa::euler_var_t{1.0, 2.0, 3.0, 4.0, 5.0};

  u(2) = ua;
  for (zisa::int_t k = 0; k < n_vars; ++k) {
    REQUIRE(u(2, k) == ua(k));
  }

  ua = u(3);
  for (zisa::int_t k = 0; k < n_vars; ++k) {
    REQUIRE(ua(k) == 0.0);
  }
}

TEST_CASE("AllVariables; load/store", "[memory]") {

  zisa::int_t n_cells = 100;
  zisa::int_t n_cvars = zisa::euler_var_t::size();
  zisa::int_t n_avars = 2;

  auto dims = zisa::AllVariablesDimensions{n_cells, n_cvars, n_avars};

  auto all_vars_store = zisa::AllVariables(dims);

  auto rd = std::random_device();
  auto gen = std::mt19937(rd());
  auto dist = std::uniform_real_distribution(-2.0, 100.0);

  for (zisa::int_t i = 0; i < all_vars_store.size(); ++i) {
    all_vars_store[i] = dist(gen);
  }

  auto filename = std::string("__unit_tests-all_variables.h5");

  { // Only one (1) HDF5 write may exist at a time.
    auto writer = zisa::HDF5SerialWriter(filename);
    zisa::save(writer, all_vars_store, zisa::all_labels<zisa::euler_var_t>());
  }

  auto reader = zisa::HDF5SerialReader(filename);
  auto all_vars_load
      = zisa::AllVariables::load(reader, zisa::all_labels<zisa::euler_var_t>());

  REQUIRE(all_vars_load.dims() == all_vars_store.dims());

  for (zisa::int_t i = 0; i < all_vars_store.size(); ++i) {
    INFO(string_format("i = %d", i));
    REQUIRE(all_vars_load[i] == all_vars_store[i]);
  }
}
