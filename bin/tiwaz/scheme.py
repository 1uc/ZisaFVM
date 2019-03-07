import json
import os

class Subsection(dict):
    def __init__(self, d):
        super().__init__(d)

    def absolute_paths(self, dirname):
        pass

class Scheme:
    def __init__(self, subsections):
        self.subsections = subsections

    def __getitem__(self, item):
        return self.subsections[item]

    def __contains__(self, item):
        return item in self.subsections

    def __iter__(self):
        return iter(self.subsections.values())

    def grid_filename(self):
        return self["grid"]["grid"]["file"]

    def absolute_paths(self, dirname):
        for subsection in self.subsections.values():
            subsection.absolute_paths(dirname)

    def save(self, filename):
        d = dict()
        for s in self.subsections.values():
            d.update(s)

        with open(filename, "w") as f:
            json.dump(d, f, indent=4)


class Experiment(Subsection):
    def __init__(self, name):
        super().__init__({
            "experiment": {
                "name": name
            }
        })

    def short_id(self):
        return self["experiment"]["name"]


class Euler(Subsection):
    def __init__(self, eos, gravity):
        super().__init__({"euler": {**eos, **gravity}})


class IdealGasEOS(Subsection):
    def __init__(self, gamma, r_gas):
        super().__init__({
            "eos": {
                "gamma": gamma, "specific-gas-constant": r_gas
            }
        })


class FluxBC(Subsection):
    def __init__(self, mode):
        super().__init__({
            "flux-bc": {
                "mode": mode
            }
        })


class PolytropeGravity(Subsection):
    def __init__(self):
        super().__init__({
            "gravity": {
                "mode": "polytrope",
                "polytrope": {
                    "rhoC": 1.0,
                    "K": 1.0,
                    "G": 1.0
                }
            }
        })


class Quadrature(Subsection):
    def __init__(self, volume_deg, edge_deg=None):
        edge_deg = edge_deg if edge_deg else volume_deg

        super().__init__({
            "quadrature": {
                "__comment" : "Specify the quadrature degree (not #points).",
                "edge": edge_deg,
                "volume": volume_deg
            }
        })


class WellBalancing(Subsection):
    def __init__(self, mode):
        super().__init__({
            "well-balancing": {
                "mode": mode
            }
        })

    def short_id(self):
        return self["well-balancing"]["mode"]


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
            "reconstruction": {
                "mode": rc,

                "orders": orders,
                "biases": biases,
                "overfit_factors": overfit_factors,
                "linear_weights": linear_weights,

                "smoothness_indicator": {
                    "epsilon": 1e-12,
                    "exponent": 4
                }
            }
        })

    def short_id(self):
        orders = self["reconstruction"]["orders"]
        return "o"+"".join(str(k) for k in orders)


class ODE(Subsection):
    def __init__(self, method, cfl_number=None):

        if not cfl_number:
            cfl_number = self.default_cfl_number(method)

        super().__init__({
            "ode": {
                "solver": method,
                "cfl_number": cfl_number
            },
        })

    def default_cfl_number(self, method):
        return 0.85

    def short_id(self):
        return self["ode"]["solver"]


class Time(Subsection):
    def __init__(self, **kwargs):
        super().__init__({"time": kwargs})


class IO(Subsection):
    def __init__(self, mode, experiment):
        super().__init__({
            "io": {
                "steps_per_frame": 1
            }
        })

        if mode == "opengl":
            self.activate_opengl()

        if mode == "hdf5":
            self.activate_hdf5(experiment)

    def activate_hdf5(self, experiment):
        self["io"]["mode"] = "hdf5"
        self["io"]["hdf5"] = {
            "filename": {
                "stem": "{:}_data-".format(experiment),
                "pattern": "%04d",
                "suffix": ".h5"
            }
        }

    def activate_opengl(self):
        self["io"]["mode"] = "opengl"
        self["io"]["mode"] = {
            "delay": "0ms"
        }

    def absolute_paths(self, dirname):
        io = self["io"]

        if "hdf5" in io:
            filename = io["hdf5"]["filename"]["stem"]
            if not os.path.isabs(filename):
                io["hdf5"]["filename"]["stem"] = os.path.join(dirname, filename)


class Grid(Subsection):
    def __init__(self, grid_name):
        super().__init__({
            "grid": {"file": grid_name}
        })

    def absolute_paths(self, dirname):
        filename = self["grid"]["file"]
        if not os.path.isabs(filename):
            self["grid"]["file"] = os.path.join(dirname, filename)
