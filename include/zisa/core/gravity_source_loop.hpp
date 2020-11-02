#ifndef GRAVITY_SOURCE_LOOP_H_F39RG
#define GRAVITY_SOURCE_LOOP_H_F39RG

#include <zisa/loops/for_each.hpp>
#include <zisa/model/all_variables.hpp>
#include <zisa/reconstruction/global_reconstruction.hpp>

namespace zisa {

template <class Equilibrium, class RC, class LEOS, class Gravity, class Scaling>
class GravitySourceLoop : public RateOfChange {
private:
  using euler_t = Euler;
  using leos_t = LEOS;
  using gravity_t = Gravity;
  using cvars_t = euler_var_t;
  using grc_t = EulerGlobalReconstruction<Equilibrium, RC, Scaling>;

public:
  GravitySourceLoop(std::shared_ptr<Grid> grid,
                    std::shared_ptr<leos_t> local_eos,
                    std::shared_ptr<gravity_t> gravity,
                    std::shared_ptr<grc_t> global_reconstruction)
      : local_eos(std::move(local_eos)),
        gravity(std::move(gravity)),
        grid(std::move(grid)),
        global_reconstruction(std::move(global_reconstruction)) {}

  virtual void compute(AllVariables &tendency,
                       const AllVariables & /* current_state */,
                       double /* t */) const override {

    auto f = [this, &tendency](int_t i, const Cell &cell) {
      const auto &rc = (*global_reconstruction)(i);
      const auto &eos = *(*local_eos)(i);

      auto x_cell = grid->cell_centers(i);

      // equilibrium terms
      auto s = cvars_t(0.0);
      for (int_t k = 0; k < grid->max_neighbours; ++k) {
        auto face = grid->face(i, k);

        auto s_eq = [&x_cell, &rc, &eos, &face](const XYZ &x) {
          auto [u_eq, w_eq] = rc.background(x);
          auto p_eq = w_eq.p;

          auto n = unit_outward_normal(face, x_cell);
          auto s = cvars_t{0.0, p_eq * n[0], p_eq * n[1], p_eq * n[2], 0.0};

          return s;
        };

        s += quadrature(face, s_eq);
      }

      // delta terms
      auto s_delta = [&rc, &gravity = *this->gravity](const XYZ &x) {
        static_assert(XYZ::size() == 3);

        auto [u_eq, _] = rc.background(x);
        auto du = rc.delta(x);
        auto u = cvars_t(u_eq + du);

        auto drho = du[0];
        auto mv = momentum(u);
        auto grad_phi = gravity.grad_phi(x);
        static_assert(decltype(grad_phi)::size() == 3);

        auto s = cvars_t{0.0,
                         -drho * grad_phi[0],
                         -drho * grad_phi[1],
                         -drho * grad_phi[2],
                         -zisa::dot(mv, grad_phi)};

        return s;
      };

      s += quadrature(cell, s_delta);
      tendency.cvars(i) += s / volume(cell);
    };

    zisa::for_each(cells(*grid), f);
  }

  virtual std::string str() const override {
    return "Gravity source term:\n" + indent_block(1, type_name<Equilibrium>())
           + "\n" + indent_block(1, type_name<RC>());
  }

private:
  std::shared_ptr<leos_t> local_eos;
  std::shared_ptr<gravity_t> gravity;
  std::shared_ptr<Grid> grid;
  std::shared_ptr<grc_t> global_reconstruction;
};

template <class RC, class EOS, class Gravity, class Scaling>
class GravitySourceLoop<NoEquilibrium, RC, EOS, Gravity, Scaling>
    : public RateOfChange {

private:
  using euler_t = Euler;
  using eos_t = EOS;
  using gravity_t = Gravity;
  using cvars_t = euler_var_t;
  using grc_t = EulerGlobalReconstruction<NoEquilibrium, RC, Scaling>;

public:
  GravitySourceLoop(std::shared_ptr<Grid> grid,
                    std::shared_ptr<eos_t>,
                    std::shared_ptr<gravity_t> gravity,
                    std::shared_ptr<grc_t> global_reconstruction)
      : gravity(std::move(gravity)),
        grid(std::move(grid)),
        global_reconstruction(std::move(global_reconstruction)) {}

  virtual void compute(AllVariables &tendency,
                       const AllVariables & /* current_state */,
                       double /* t */) const override {

    auto f = [this, &tendency](int_t i, const Cell &cell) {
      const auto &rc = (*global_reconstruction)(i);

      auto s = [&rc, &gravity = *this->gravity](const XYZ &x) {
        static_assert(XYZ::size() == 3);

        auto u = rc(x);
        auto rho = u[0];
        auto mv = momentum(u);
        auto grad_phi = gravity.grad_phi(x);
        static_assert(decltype(grad_phi)::size() == 3);

        return cvars_t{0.0,
                       -rho * grad_phi[0],
                       -rho * grad_phi[1],
                       -rho * grad_phi[2],
                       -zisa::dot(mv, grad_phi)};
      };

      tendency.cvars(i) += average(cell, s);
    };

    zisa::for_each(cells(*grid), f);
  }

  virtual std::string str() const override {
    return "Gravity source term:\n"
           + indent_block(1, type_name<NoEquilibrium>()) + "\n"
           + indent_block(1, type_name<RC>());
  }

private:
  std::shared_ptr<gravity_t> gravity;
  std::shared_ptr<Grid> grid;
  std::shared_ptr<grc_t> global_reconstruction;
};

} // namespace zisa
#endif /* end of include guard */
