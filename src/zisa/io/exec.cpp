/* C++11 version of 'exec'.
 *
 * Authors: Luc Grosheintz <forbugrep@zoho.com>
 *    Date: 2016-01-14
 */

#include <cstdlib>
#include <iostream>
#include <memory>

#include <zisa/io/exec.hpp>

namespace zisa {

std::string exec(const std::string &cmd) {
  std::shared_ptr<FILE> pipe(popen(cmd.c_str(), "r"), pclose);
  if (!pipe) {
    std::cerr << "Failed to open a pipe.\n";
    exit(EXIT_FAILURE);
  }

  std::string result = "";

  // Read the buffer in 128 byte chunks.
  char read_buffer[128];
  while (!feof(pipe.get())) {
    if (fgets(read_buffer, 128, pipe.get()) != NULL) {
      result += read_buffer;
    }
  }

  return result;
}

} // namespace zisa
