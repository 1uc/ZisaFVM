#ifndef SANITY_CHECK_FOR_H_UX4Z8
#define SANITY_CHECK_FOR_H_UX4Z8

#include <zisa/config.hpp>
#include <zisa/model/all_variables.hpp>
#include <zisa/model/sanity_check.hpp>

namespace zisa {

template <class Model>
class SanityCheckFor : public SanityCheck {
public:
  SanityCheckFor(const Model &model) : model(model) {}

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
  Model model;
};

} // namespace zisa

#endif /* end of include guard */
