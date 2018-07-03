/* Input parameters.
 *
 */
#ifndef INPUT_PARAMETERS_H_L08MP
#define INPUT_PARAMETERS_H_L08MP

#include <string>
#include <nlohmann/json.hpp>

namespace zisa {

class InputParameters {
private:
  nlohmann::json json_;

public:
  InputParameters(const std::string &filename);

  auto operator[](const std::string &key) const -> decltype(json_["key"]);
  auto operator[](const std::string &key) -> decltype(json_["key"]);

};

} // namespace zisa

#endif /* end of include guard */
