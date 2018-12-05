/* Simple progress bar.
 *
 * Authors: Luc Grosheintz <forbugrep@zoho.com>
 *    Date: 2015-01-12
 */
#include <zisa/datetime.hpp>
#include <zisa/io/progress_bar.hpp>

namespace zisa {
ProgressBar::ProgressBar() { t_start = zisa::current_time_stamp(); }

ProgressBar::ProgressBar(int unit_size) : unit_size(unit_size) {
  t_start = zisa::current_time_stamp();
}

void ProgressBar::reset(void) {
  t_start = zisa::current_time_stamp();
  step = 0;
}
} // namespace zisa
