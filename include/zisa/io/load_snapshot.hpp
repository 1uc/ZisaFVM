// SPDX-License-Identifier: MIT
// Copyright (c) 2021 ETH Zurich, Luc Grosheintz-Laval

#ifndef ZISA_LOAD_SNAPSHOT_HPP_BCZZLK
#define ZISA_LOAD_SNAPSHOT_HPP_BCZZLK

#include <zisa/config.hpp>
#include <zisa/io/data_source.hpp>
#include <zisa/io/file_name_generator.hpp>
#include <zisa/io/hdf5_serial_writer.hpp>
#include <zisa/model/all_variables_fwd.hpp>
#include <zisa/model/load_full_state.hpp>
#include <zisa/model/models.hpp>
#include <zisa/ode/simulation_clock.hpp>

namespace zisa {

class LoadSnapshot : public DataSource {
public:
  LoadSnapshot(std::shared_ptr<FNG> fng);

protected:
  virtual void do_load(AllVariables &all_variables,
                       SimulationClock &simulation_clock) override;

protected:
  virtual std::unique_ptr<HDF5Reader> pick_reader(const std::string &file_name)
      = 0;

private:
  std::shared_ptr<FNG> fng;
};

class SerialLoadSnapshot : public LoadSnapshot {
private:
  using super = LoadSnapshot;

public:
  using super::super;

protected:
  virtual std::unique_ptr<HDF5Reader>
  pick_reader(const std::string &file_name) override;
};

}

#endif // ZISA_LOAD_SNAPSHOT_HPP
