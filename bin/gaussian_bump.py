#! /usr/bin/env python3

import os
import shutil
import glob

import tiwaz
import tiwaz.scheme as sc

from tiwaz.launch_params import all_combinations, pointwise_combinations
from tiwaz.launch_params import build_zisa

from tiwaz.post_process import load_results

from tiwaz.cli_parser import default_cli_parser
from tiwaz.launch_job import launch_all
from tiwaz.latex_tables import write_convergence_table
from tiwaz.scatter_plot import plot_visual_convergence

class GaussianBumpExperiment(sc.Subsection):
    def __init__(self, amplitude, width):
        super().__init__({
            "name": "gaussian_bump",

            "initial_conditions": {
                "amplitude": amplitude,
                "width": width
            }
        })

    def short_id(self):
        amp = self["initial_conditions"]["amplitude"]
        return self["name"] + "_amp{:.2e}".format(amp)


amplitudes = [0.0]
width = 0.05

eos = sc.IdealGasEOS(gamma=2.0, r_gas=1.0)
gravity = sc.PolytropeGravity()
euler = sc.Euler(eos, gravity)

time = sc.Time(t_end=0.09)
io = sc.IO("hdf5", "gaussian_bump", n_snapshots=1)

def grid_name(level):
    return "grids/polytrope-{:}.msh".format(level)

coarse_grid_levels = list(range(0, 4))
coarse_grid_names = [grid_name(level) for level in coarse_grid_levels]

coarse_grid_choices = {
    "grid": [
        sc.Grid(grid_name(l), l) for l in coarse_grid_levels
    ]
}
coarse_grids = all_combinations(coarse_grid_choices)
reference_grid = sc.Grid("grids/polytrope-4.msh", 4)

independent_choices = {
    "euler": [euler],

    "flux-bc": [sc.FluxBC("isentropic")],

    "well-balancing": [sc.WellBalancing("constant"),
                       sc.WellBalancing("isentropic")],

    "io": [io],

    "time": [time],
}

# dependent_choices = {
#     "reconstruction": [
#         # sc.Reconstruction("CWENO-AO", [1]),
#         sc.Reconstruction("CWENO-AO", [2, 2, 2, 2])
#         #sc.Reconstruction("CWENO-AO", [3, 2, 2, 2])
#     ],
#
#     "ode": [
#         # sc.ODE("ForwardEuler"),
#         sc.ODE("SSP3")
#         # sc.ODE("SSP3")
#     ],
#
#     "quadrature": [
#         # sc.Quadrature(1),
#         sc.Quadrature(1)
#         # sc.Quadrature(3)
#     ]
# }

dependent_choices = {
    "reconstruction": [
        sc.Reconstruction("CWENO-AO", [1]),
        sc.Reconstruction("CWENO-AO", [2, 2, 2, 2], overfit_factors=[3.0, 2.0, 2.0, 2.0]),
        sc.Reconstruction("CWENO-AO", [3, 2, 2, 2]),
        sc.Reconstruction("CWENO-AO", [4, 2, 2, 2])
    ],

    "ode": [
        sc.ODE("ForwardEuler"),
        sc.ODE("SSP2"),
        sc.ODE("SSP3"),
        sc.ODE("SSP3")
    ],

    "quadrature": [
        sc.Quadrature(1),
        sc.Quadrature(2),
        sc.Quadrature(3),
        sc.Quadrature(4)
    ]
}

reference_choices = {
    "euler": [euler],
    "flux-bc": [sc.FluxBC("isentropic")],
    "well-balancing": [sc.WellBalancing("isentropic")],
    "time": [time],
    "io": [io],
    "reconstruction": [sc.Reconstruction("CWENO-AO", [5, 2, 2, 2])],
    "ode": [sc.ODE("SSP3")],
    "quadrature": [sc.Quadrature(4)],
    "grid": [reference_grid],
    "reference": [sc.Reference("isentropic", coarse_grid_names)]
}

base_choices = all_combinations(independent_choices)
model_choices = base_choices.product(pointwise_combinations(dependent_choices))

coarse_runs_ = model_choices.product(coarse_grids)
reference_runs_ = all_combinations(reference_choices)

def make_runs(amplitude):
    coarse_runs = coarse_runs_.product([{
        "experiment": GaussianBumpExperiment(amplitude, width)
    }])

    reference_runs = reference_runs_.product([{
        "experiment": GaussianBumpExperiment(amplitude, width)
    }])

    coarse_runs = [sc.Scheme(choice) for choice in coarse_runs]
    reference_runs = [sc.Scheme(choice) for choice in reference_runs]

    return coarse_runs, reference_runs

all_runs = [make_runs(amp) for amp in amplitudes]

def post_process(coarse_runs, reference_run):
    results, columns = load_results(coarse_runs, reference_run)
    labels = TableLabels()

    filename = coarse_runs[0]["experiment"].short_id()
    write_convergence_table(results, columns, labels, filename)
    plot_visual_convergence(results, columns, labels, filename)

class TableLabels:
    def __call__(self, col):
        return " ".join(str(v) for v in col.values())

def main():
    parser = default_cli_parser("'gaussian_bump' numerical experiment.")
    args = parser.parse_args()

    if args.run:
        build_zisa()

        for c, r in all_runs:
            launch_all(c, force=args.force)

            if args.reference:
                launch_all(r, force=args.force)

    if args.post_process:
        for c, r in all_runs:
            post_process(c, r[0])

    if args.copy_to_paper:
        dir = "${HOME}/git/papers/LucGrosheintz/papers/unstructured_well_balancing/img/gaussian_bump"
        dir = os.path.expandvars(dir)

        for c, _ in all_runs:
            stem = c[0]["experiment"].short_id()
            files = glob.glob(stem + "*.tex")

            for f in files:
                shutil.copy(f, os.path.join(dir, os.path.basename(f)))

if __name__ == "__main__":
    main()







