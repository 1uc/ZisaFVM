/* Simple Linux backtrace with line numbers.
 *
 * Authors: Luc Grosheintz <forbugrep@zoho.com>
 *    Date: 2016-01-14
 */
#if defined(__GNUC__) || defined(__GNUG__)
#include <execinfo.h>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>

#include <zisa/io/backtrace.hpp>
#include <zisa/io/exec.hpp>
#include <zisa/utils/string_format.hpp>

namespace zisa {
std::string backtrace(void) {
  const int trace_depth = 16;

  void *trace[trace_depth];
  char **messages = nullptr;

  int trace_size = ::backtrace(trace, trace_depth);
  messages = backtrace_symbols(trace, trace_size);

  std::stringstream ss;
  ss << "[TRACE] Backtrace: \n";

  // Skip i == 0, which is this function.
  for (int i = 1; i < trace_size; ++i) {
    // The first bracket or space should be where the file name ends.
    int p = 0;
    while (messages[i][p] != '(' && messages[i][p] != ' '
           && messages[i][p] != 0) {
      ++p;
    }

    char command_to_execute[256];
    snprintf(command_to_execute,
             256,
             "addr2line %p -e %.*s",
             trace[i],
             p,
             messages[i]);

    ss << " ||  " << zisa::exec(command_to_execute);
  }
  ss << "[=====] \n";

  return ss.str();
}
} // namespace zisa
#else  // any compiler other than GNU C++
namespace zisa {
std::string backtrace(void) {
  return "No backtrace available on this platform!";
}
} // namespace zisa
#endif // GNU C++
