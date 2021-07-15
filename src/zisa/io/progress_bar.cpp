// SPDX-License-Identifier: MIT
// Copyright (c) 2021 ETH Zurich, Luc Grosheintz-Laval

#include <zisa/datetime.hpp>
#include <zisa/io/progress_bar.hpp>

namespace zisa {
SerialProgressBar::SerialProgressBar() { t_start = zisa::current_time_stamp(); }

SerialProgressBar::SerialProgressBar(int unit_size) : unit_size(unit_size) {
  t_start = zisa::current_time_stamp();
}

void SerialProgressBar::write_progress(std::ostream &out,
                                       const std::string &msg) {
  // if not active, do nothing.
  if (!is_active) {
    return;
  }

  ++step;

  if (step % unit_size == 0) {
    out << unit_symbol;
  }
  if (step % (5 * unit_size) == 0) {
    out << " ";
  }
  if (step % (10 * unit_size) == 0) {
    auto elapsed = elapsed_seconds_since(t_start);
    auto time_formatted = duration_format(elapsed);

    out << " " << msg;
    out << " (";
    out << time_formatted;
    out << "; " << step << ")\n";
  }
  out.flush();
}

void SerialProgressBar::activate() { is_active = true; }
void SerialProgressBar::deactivate() { is_active = false; }

void SerialProgressBar::reset() {
  t_start = zisa::current_time_stamp();
  step = 0;
}

} // namespace zisa
