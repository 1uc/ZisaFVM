// SPDX-License-Identifier: MIT
// Copyright (c) 2021 ETH Zurich, Luc Grosheintz-Laval

#ifndef ZISA_HEATING_HPP
#define ZISA_HEATING_HPP

#include <zisa/config.hpp>

#include <functional>
#include <memory>

#include <zisa/cli/input_parameters.hpp>
#include <zisa/grid/grid.hpp>
#include <zisa/math/quadrature.hpp>

namespace zisa {

template <class GRC>
class Heating : public RateOfChange {
public:
  Heating(std::shared_ptr<Grid> grid,
          std::shared_ptr<GRC> grc,
          std::function<double(const XYZ &x)> heating_rate)
      : grid(std::move(grid)),
        grc(std::move(grc)),
        heating_rate(std::move(heating_rate)) {}

  virtual void compute(AllVariables &tendency,
                       const AllVariables & /* current_state */,
                       double /* t */) const override {

    auto &dudt = tendency.cvars;

    zisa::for_each(cells(*grid), [this, &dudt](int_t i, const Cell &cell) {
      auto f = [this, i](const XYZ &x) {
        auto rho
            = [this, i](const XYZ &x) -> double { return (*grc)(i)(x)[0]; };

        return rho(x) * heating_rate(x);
      };

      dudt(i, 4) += average(cell, f);
    });
  }

  virtual std::string str() const override { return "Constant heating rate."; }

private:
  std::shared_ptr<Grid> grid;
  std::shared_ptr<GRC> grc;
  std::function<double(const XYZ &x)> heating_rate;
};

template <class GRC>
std::shared_ptr<RateOfChange>
make_heating_source(const std::shared_ptr<Grid> &grid,
                    const std::shared_ptr<GRC> &grc,
                    const InputParameters &params) {

  if (has_key(params, "heating")) {
    double epsilon = params["heating"]["rate"];
    double r0 = params["heating"]["lower_boundary"];
    double r1 = params["heating"]["upper_boundary"];

    auto xi = [r0, r1](const XYZ &x) {
      double r = zisa::norm(x);
      return (r0 <= r) && (r <= r1);
    };

    auto heating_rate
        = [grid, epsilon, xi](const XYZ &x) { return epsilon * xi(x); };

    return std::make_shared<Heating<GRC>>(grid, grc, heating_rate);
  }

  return nullptr;
}

}

#endif // ZISA_HEATING_HPP
