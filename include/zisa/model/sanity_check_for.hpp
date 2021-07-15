// SPDX-License-Identifier: MIT
// Copyright (c) 2021 ETH Zurich, Luc Grosheintz-Laval

#ifndef SANITY_CHECK_FOR_H_UX4Z8
#define SANITY_CHECK_FOR_H_UX4Z8

#include <zisa/config.hpp>

#include <memory>
#include <zisa/model/all_variables.hpp>
#include <zisa/model/sanity_check.hpp>

namespace zisa {

template <class Model>
class SanityCheckFor : public SanityCheck {
private:
  using cvars_t = typename Model::cvars_t;

public:
  explicit SanityCheckFor(std::shared_ptr<Model> model)
      : model(std::move(model)) {}

  bool operator()(const AllVariables &all_variables) const override {
    auto &cvars = all_variables.cvars;

    for (int_t i = 0; i < cvars.shape(0); ++i) {
      auto u = cvars_t(cvars(i));

      if (!isplausible(u)) {
#if ZISA_HAS_MPI == 1
        std::cout << string_format("[PE %d] (%d) %s\n",
                                   zisa::mpi::rank(MPI_COMM_WORLD),
                                   i,
                                   format_as_list(u).c_str());
#else
        std::cout << string_format("(%d) [", i) << u << "]\n";
#endif
        return false;
      }
    }

    return true;
  }

private:
  std::shared_ptr<Model> model;
};

} // namespace zisa

#endif /* end of include guard */
