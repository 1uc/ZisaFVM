#ifndef SANITY_CHECK_FOR_H_UX4Z8
#define SANITY_CHECK_FOR_H_UX4Z8

#include <zisa/config.hpp>

#include <memory>
#include <zisa/model/all_variables.hpp>
#include <zisa/model/sanity_check.hpp>

namespace zisa {

template <class Model>
class SanityCheckFor : public SanityCheck {
public:
  explicit SanityCheckFor(std::shared_ptr<Model> model)
      : model(std::move(model)) {}

  bool operator()(const AllVariables &all_variables) const override {
    using cvars_t = typename Model::cvars_t;
    auto &cvars = all_variables.cvars;

    for (int_t i = 0; i < cvars.shape(0); ++i) {
      auto u = cvars_t(cvars(i));

      if (!isplausible(u)) {
        std::cout << string_format("(%d) [", i) << u << "]\n";
        return false;
      }
    }

    return true;
  }

private:
  std::shared_ptr<Model> model;
};

template <>
class SanityCheckFor<Euler<JankaEOS, RadialGravity>> : public SanityCheck {
private:
  using euler_t = Euler<JankaEOS, RadialGravity>;

public:
  explicit SanityCheckFor(std::shared_ptr<euler_t> euler)
      : euler(std::move(euler)) {}

  inline bool operator()(const AllVariables &all_variables) const override {
    auto &cvars = all_variables.cvars;

    for (int_t i = 0; i < cvars.shape(0); ++i) {
      auto u = euler_var_t(cvars(i));

      bool is_rho_positive = (u[0] > 0.0);
      //      bool is_super_polytropic
      //          = (u[4] >= (1.0 - 1e-9) * euler->eos.polytropic_energy(u[0]));
      bool is_super_polytropic = true; // turn off for now.
      bool is_real = zisa::isreal(u);

      if (!(is_rho_positive && is_super_polytropic && is_real)) {
        std::cout << string_format("(%d) [", i) << u << "]\n";
        return false;
      }
    }

    return true;
  }

private:
  std::shared_ptr<euler_t> euler;
};

} // namespace zisa

#endif /* end of include guard */
