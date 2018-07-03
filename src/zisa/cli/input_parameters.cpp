/* Input parameters.
 *
 */

#include <fstream>

#include <zisa/utils/logging.hpp>
#include <zisa/cli/input_parameters.hpp>
#include <zisa/utils/string_format.hpp>

namespace zisa {

InputParameters::InputParameters(const std::string &filename)
{
  auto is = std::ifstream(filename);
  LOG_ERR_IF(!is, string_format("Could not read file. [%s]", filename.c_str()));

  is >> json_;
}

auto InputParameters::operator[](const std::string &key) const
    -> decltype(json_["key"]) {
  return json_[key];
}

auto InputParameters::operator[](const std::string &key)
    -> decltype(json_["key"]) {
  return json_[key];
}

} // namespace zisa
