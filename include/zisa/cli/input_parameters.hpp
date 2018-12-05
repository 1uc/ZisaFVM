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

  inline decltype(auto) operator[](const std::string &key) const { return json_[key]; }
  inline decltype(auto) operator[](const std::string &key) { return json_[key]; }

};

} // namespace zisa

#endif /* end of include guard */
