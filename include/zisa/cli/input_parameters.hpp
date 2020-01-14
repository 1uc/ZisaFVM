/* Input parameters.
 *
 */
#ifndef INPUT_PARAMETERS_H_L08MP
#define INPUT_PARAMETERS_H_L08MP

#include <nlohmann/json.hpp>
#include <string>

#include <zisa/utils/has_key.hpp>

namespace zisa {

class InputParameters {
private:
  nlohmann::json json_;

public:
  InputParameters(const std::string &filename);

  inline decltype(auto) operator[](const std::string &key) const {
    return json_[key];
  }
  inline decltype(auto) operator[](const std::string &key) {
    return json_[key];
  }

  bool has_key(const std::string &key) const;
};

inline bool has_key(const InputParameters &params, const std::string &key) {
  return params.has_key(key);
}

inline bool is_mpi(const InputParameters &params) {
  return has_key(params, "parallelization")
         && has_key(params["parallelization"], "mode")
         && params["parallelization"]["mode"] == "mpi";
}

} // namespace zisa

#endif /* end of include guard */
