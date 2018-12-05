/* Input parameters.
 *
 */

#include <fstream>

#include <zisa/cli/input_parameters.hpp>
#include <zisa/utils/logging.hpp>
#include <zisa/utils/string_format.hpp>

namespace zisa {

InputParameters::InputParameters(const std::string &filename) {
  auto is = std::ifstream(filename);
  LOG_ERR_IF(!is, string_format("Could not read file. [%s]", filename.c_str()));

  is >> json_;
}

bool InputParameters::has_key(const std::string &key) const {
  return json_.find(key) != json_.end();
}

} // namespace zisa
