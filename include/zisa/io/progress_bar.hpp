/* Simple CLI progress bar.
 *
 * Authors: Luc Grosheintz <forbugrep@zoho.com>
 *    Date: 2015-01-12
 */
#ifndef PROGRESS_BAR_H_JSW0ELA2
#define PROGRESS_BAR_H_JSW0ELA2
#include <iostream>

#include <zisa/config.hpp>
#include <zisa/datetime.hpp>

namespace zisa {

class ProgressBar {
public:
  ProgressBar();
  ProgressBar(int unit_size);

  /// Write progress into stream.
  template <class Stream>
  void write_progress(Stream &out, const std::string &msg) {
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

  inline void activate() { is_active = true; }
  inline void deactivate() { is_active = false; }
  void reset(void);

private:
  int unit_size = 100;
  char unit_symbol = '.';
  zisa::time_stamp_t t_start;
  bool is_active = true;
  int step = 0;
};

} // namespace zisa

#endif /* end of include guard: PROGRESS_BAR_H_JSW0ELA2 */
