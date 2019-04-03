import json

class Subsection(dict):
    def __init__(self, d):
        super().__init__(d)

class Scheme:
    def __init__(self, subsections):
        self.subsections = subsections

    def __getitem__(self, item):
        return self.subsections[item]

    def __contains__(self, item):
        return item in self.subsections

    def __iter__(self):
        return iter(self.subsections.values())

    def order(self):
        return self["reconstruction"]["orders"][0]

    def well_balancing(self):
        return self["well-balancing"]["mode"]

    def grid_level(self):
        return self["grid"]["level"]

    def grid_filename(self):
        return self["grid"]["file"]

    def folder_name(self):
        subsections = ['experiment',
                       'reconstruction',
                       'order', 'ode',
                       'well-balancing', 'grid']

        subsections = [s for s in subsections if s in self]

        name = "_".join([self[key].short_id() for key in subsections])

        return name


    def save(self, filename):
        with open(filename, "w") as f:
            json.dump(self.subsections, f, indent=4)


class Euler(Subsection):
    def __init__(self, eos, gravity):
        super().__init__({"eos": eos, "gravity": gravity})


class IdealGasEOS(Subsection):
    def __init__(self, gamma, r_gas):
        super().__init__({
            "gamma": gamma, "specific-gas-constant": r_gas
        })


class FluxBC(Subsection):
    def __init__(self, mode):
        super().__init__({
            "mode": mode
        })


class PolytropeGravity(Subsection):
    def __init__(self, rhoC=1.0, K=1.0, G=1.0):
        super().__init__({
            "mode": "polytrope",
            "rhoC": rhoC,
            "K": K,
            "G": G
        })

class PolytropeGravityWithJump(Subsection):
    def __init__(self, r_crit=0.25, rhoC=1.0, K_inner=1.0, K_outer=1.0, G=1.0):
        super().__init__({
            "mode": "polytrope_with_jump",
            "r_crit": r_crit,
            "rhoC": rhoC,
            "K_inner": K_inner,
            "K_outer": K_outer,
            "G": G
        })


class Quadrature(Subsection):
    def __init__(self, volume_deg, edge_deg=None):
        edge_deg = edge_deg if edge_deg else volume_deg

        super().__init__({
            "__comment" : "Specify the quadrature degree (not #points).",
            "edge": edge_deg,
            "volume": volume_deg
        })


class WellBalancing(Subsection):
    def __init__(self, mode):
        super().__init__({
            "mode": mode
        })

    def short_id(self):
        return self["mode"]


class Reconstruction(Subsection):
    def __init__(self, rc, orders, biases=None, overfit_factors=None, linear_weights=None):
        n_stencils = len(orders)

        if not biases:
            biases = ["c" if i == 0 else "b" for i in range(n_stencils)]

        if not overfit_factors:
            overfit_factors = [2.0 if i == 0 else 1.5 for i in range(n_stencils)]

        if not linear_weights:
            linear_weights = [100 if i == 0 else 1 for i in range(n_stencils)]

        super().__init__({
            "mode": rc,

            "orders": orders,
            "biases": biases,
            "overfit_factors": overfit_factors,
            "linear_weights": linear_weights,

            "smoothness_indicator": {
                "epsilon": 1e-10,
                "exponent": 4
            }
        })

    def short_id(self):
        orders = self["orders"]
        return "o"+"".join(str(k) for k in orders)


class ODE(Subsection):
    def __init__(self, method, cfl_number=None):

        if not cfl_number:
            cfl_number = self.default_cfl_number(method)

        super().__init__({
            "solver": method,
            "cfl_number": cfl_number
        })

    def default_cfl_number(self, method):
        return 0.85

    def short_id(self):
        return self["solver"]


class Time(Subsection):
    def __init__(self, **kwargs):
        super().__init__(kwargs)


class IO(Subsection):
    def __init__(self, mode, experiment, **kwargs):
        super().__init__(kwargs)

        self["filename"] = {
            "stem": "{:}_data-".format(experiment),
            "pattern": "%04d",
            "suffix": ".h5"
        }

        if mode == "opengl":
            self.activate_opengl()

        if mode == "hdf5":
            self.activate_hdf5(experiment)

    def activate_hdf5(self, experiment):
        self["mode"] = "hdf5"

    def activate_opengl(self):
        self["mode"] = "opengl"
        self["opengl"] = {
            "delay": "0ms"
        }


class Grid(Subsection):
    def __init__(self, grid_name, level):
        super().__init__({
            "file": grid_name,
             "level": level
        })

    def short_id(self):
        l = self["level"]
        return "L" + str(l)


class Reference(Subsection):
    def __init__(self, equilibrium, coarse_grid_names):
        super().__init__({
            "equilibrium": equilibrium,
            "coarse_grids": coarse_grid_names
        })