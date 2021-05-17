#ifndef PROGRESS_BAR_H_JSW0ELA2
#define PROGRESS_BAR_H_JSW0ELA2
#include <iostream>

#include <zisa/config.hpp>
#include <zisa/datetime.hpp>

namespace zisa {

class ProgressBar {
public:
  virtual ~ProgressBar() = default;

  virtual void write_progress(std::ostream &out, const std::string &msg) = 0;
  virtual void activate() = 0;
  virtual void deactivate() = 0;
  virtual void reset() = 0;
};

class SerialProgressBar : public ProgressBar {
public:
  SerialProgressBar();
  explicit SerialProgressBar(int unit_size);

  /// Write progress into stream.
  void write_progress(std::ostream &out, const std::string &msg) override;

  void activate() override;
  void deactivate() override;
  void reset() override;

private:
  int unit_size = 100;
  char unit_symbol = '.';
  zisa::time_stamp_t t_start;
  bool is_active = true;
  int step = 0;
};

} // namespace zisa

#endif /* end of include guard: PROGRESS_BAR_H_JSW0ELA2 */
