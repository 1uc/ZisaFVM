#!/usr/bin/env python3

import os
import shutil
import glob

import numpy as np
import matplotlib.pyplot as plt
import h5py

import tiwaz
import tiwaz.scheme as sc
import tiwaz.gmsh as gmsh

from tiwaz.launch_params import all_combinations, pointwise_combinations
from tiwaz.launch_params import build_zisa

from tiwaz.cli_parser import default_cli_parser
from tiwaz.launch_job import launch_all, restart_all
from tiwaz.scatter_plot import ScatterPlot, plot_visual_convergence
from tiwaz.launch_params import folder_name
from tiwaz.post_process import load_data, load_grid
from tiwaz.post_process import find_data_files, find_last_data_file
from tiwaz.post_process import find_steady_state_file
from tiwaz.utils import read_txt, write_txt
from tiwaz.tri_plot import TriPlot
from tiwaz.hd2d import hd2d

polytropic_index_n = 3.0
polytropic_gamma = (polytropic_index_n + 1.0) / polytropic_index_n
polytrope_radius = 1.45e8
rho_center = 1e10
K_center = 4.897e14
gamma1 = 1.325
# gamma1 = 1.25

mesh_levels = list(range(0, 2))
coarse_grid_levels = list(range(1, 2))


janka_params = dict(
    rho_bounce=2e14,
    gamma1=gamma1,
    gamma2=2.5,
    gamma_thermal=1.5,
    E1=K_center / (gamma1 - 1.0),
)


class JankaExperiment(sc.Subsection):
    def __init__(self):
        super().__init__(
            {"name": "janka", "initial_conditions": {"gamma": polytropic_gamma}}
        )

    def short_id(self):
        return self["name"]


# This is just in case we want to run the same experiment with many parameters.
experiment_params = [None]

eos = sc.JankaEOS(**janka_params)
gravity = sc.GeneralPolytropeGravity(
    rho_center=rho_center,
    polytropic_index_n=polytropic_index_n,
    radius=polytrope_radius,
    G=6.674e-8,
)
euler = sc.Euler(eos, gravity)

time = sc.Time(t_end=0.5)
io = sc.IO("hdf5", "janka", steps_per_frame=100)
# io = sc.IO("opengl", "janka", steps_per_frame=1)


def grid_name_stem(l):
    return "grids/janka-{}".format(l)


def grid_name(l):
    return grid_name_stem(l) + ".msh"


def geo_grid_name(l):
    return grid_name_stem(l) + ".geo"


mesh_sizes = [(1e-3 * 0.5 ** l, 0.3e-1 * 0.5 ** l) for l in mesh_levels]
coarse_grid_names = [grid_name(level) for level in coarse_grid_levels]

coarse_grid_choices = {"grid": [sc.Grid(grid_name(l), l) for l in coarse_grid_levels]}
coarse_grids = all_combinations(coarse_grid_choices)
reference_grid = sc.Grid(grid_name(mesh_levels[-1]), mesh_levels[-1])

independent_choices = {
    "euler": [euler],
    # "flux-bc": [sc.FluxBC("isentropic")],
    "flux-bc": [sc.FluxBC("constant")],
    # "well-balancing": [ sc.WellBalancing("isentropic") ],
    "well-balancing": [sc.WellBalancing("constant")],
    # "well-balancing": [ sc.WellBalancing("constant"),
    #                     sc.WellBalancing("isentropic")],
    "io": [io],
    "time": [time],
}

dependent_choices = {
    "reconstruction": [
        # sc.Reconstruction("CWENO-AO", [1]),
        # sc.Reconstruction("CWENO-AO", [2, 2, 2, 2])
        sc.Reconstruction(
            "CWENO-AO", [3, 2, 2, 2], overfit_factors=[3.0, 2.0, 2.0, 2.0]
        )
        # sc.Reconstruction("CWENO-AO", [5, 2, 2, 2], overfit_factors = [3.0, 2.0, 2.0, 2.0])
    ],
    "ode": [
        # sc.ODE("ForwardEuler")
        sc.ODE("SSP3", cfl_number=0.8)
        # sc.ODE("SSP3")
        # sc.ODE("SSP3")
    ],
    "quadrature": [
        # sc.Quadrature(1),
        # sc.Quadrature(1)
        sc.Quadrature(1)
        # sc.Quadrature(4)
    ],
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
    "reference": [sc.Reference("isentropic", coarse_grid_names)],
}

base_choices = all_combinations(independent_choices)
model_choices = base_choices.product(pointwise_combinations(dependent_choices))

coarse_runs_ = model_choices.product(coarse_grids)
reference_runs_ = all_combinations(reference_choices)


def make_runs(params):
    coarse_runs = coarse_runs_.product([{"experiment": JankaExperiment()}])

    reference_runs = reference_runs_.product([{"experiment": JankaExperiment()}])

    coarse_runs = [sc.Scheme(choice) for choice in coarse_runs]
    reference_runs = [sc.Scheme(choice) for choice in reference_runs]

    return coarse_runs, reference_runs


all_runs = [make_runs(param) for param in experiment_params]


def post_process(coarse_runs, reference_run):
    ref_dir = os.path.expandvars(
        "${SCRATCH}/1d_reference/janka/hd2d_collapse_4Luc/data"
    )

    # a = hd2d(ref_dir, "test1d", 201)
    #
    # rho_ref = a.v[0].reshape(-1)
    # phi_ref = a.v[23].reshape(-1)
    # r_ref = a.x1.reshape(-1)
    #
    # with h5py.File("recomputed_gravity.h5", "r") as h5:
    #     phi = np.array(h5["model/gravity/phi"])[:-1]
    #     r = np.array(h5["model/gravity/radii"])[:-1]
    #     rho = np.array(h5["rho"])
    #
    # with h5py.File("grid.h5", "r") as h5:
    #     x = np.linalg.norm(np.array(h5["cell_centers"]), axis=1)
    #
    # plt.plot(r_ref, phi_ref - phi_ref[0], label="ref")
    # plt.plot(r, phi - phi[0], label="approx")
    # plt.legend()
    # plt.show()
    #
    # def avg(x):
    #     return 0.5*(x[1:] + x[:-1])
    #
    # def diff(x):
    #     return x[1:] - x[:-1]
    #
    # plt.loglog(avg(r_ref), diff(phi_ref) / diff(r_ref), label="ref")
    # plt.loglog(avg(r), diff(phi) / diff(r), label="approx")
    # plt.legend()
    # plt.show()
    #
    # plt.loglog(r_ref, rho_ref, '+', label="ref")
    # plt.loglog(x, rho, '.', label="approx")
    # plt.legend()
    # plt.show()

    for coarse_run in coarse_runs:
        coarse_dir = folder_name(coarse_run)
        coarse_grid = load_grid(coarse_dir)

        output_dir = "{}/img".format(coarse_dir)
        os.makedirs(output_dir, exist_ok=True)

        x, y = coarse_grid.cell_centers[:, 0], coarse_grid.cell_centers[:, 1]

        # u_ideal = load_data(coarse_dir + "/ic_ideal_gas_eos.h5", None)
        # u_janka = load_data(coarse_dir + "/ic_janka_eos.h5", None)
        #
        # plot = ScatterPlot()
        # plot(coarse_grid, u_ideal.cvars["rho"])
        # plot(coarse_grid, u_janka.cvars["rho"])
        # plt.title("rho")
        #
        # plot = ScatterPlot()
        # plot(coarse_grid, u_ideal.xvars["p"])
        # plot(coarse_grid, u_janka.xvars["p"])
        # plt.legend(["ideal", "janka"])
        # plt.title("p")
        # plt.show()

        data_files = find_data_files(coarse_dir)

        u0 = load_data(data_files[0], find_steady_state_file(coarse_dir))
        for k, data_file in enumerate(data_files[:]):
            u_coarse = load_data(data_file, find_steady_state_file(coarse_dir))

            t = u_coarse.time

            a = hd2d(ref_dir, "test1d", int(round(t * 1000)))
            print(t)
            print(a.time)

            rho = u_coarse.cvars["rho"]
            rho0 = u0.cvars["rho"]

            drho = u_coarse.dvars["rho"]
            vx = u_coarse.cvars["mv1"] / rho
            vy = u_coarse.cvars["mv2"] / rho
            v = (vx * x + vy * y) / (x ** 2 + y ** 2) ** 0.5
            p = u_coarse.xvars["p"]
            p0 = u0.xvars["p"]
            cs = u_coarse.xvars["cs"]

            E = u_coarse.cvars["E"]
            E_th = u_coarse.xvars["E_th"]

            # trip = TriPlot()
            # trip.color_plot(coarse_grid, E_th)
            # # trip.quiver(coarse_grid, vx, vy)
            # plt.savefig("{}/img/Eth-colorp-{:04d}.png".format(coarse_dir, k))

            plot = ScatterPlot()
            plot(coarse_grid, E_th / E)
            plot.annotate_time(u_coarse.time)
            plt.ylim([-0.001, 0.001])
            plt.savefig("{}/img/Eth-{:04d}.png".format(coarse_dir, k))

            plot = ScatterPlot()
            plot(coarse_grid, cs)
            plt.plot(a.x1, a.v[16])
            plot.annotate_time(u_coarse.time)
            plt.savefig("{}/img/cs-{:04d}.png".format(coarse_dir, k))

            plot = ScatterPlot()
            plot(coarse_grid, p)
            plt.plot(a.x1, a.v[9])
            plot.annotate_time(u_coarse.time)
            plt.savefig("{}/img/p-{:04d}.png".format(coarse_dir, k))

            plot = ScatterPlot()
            plot(coarse_grid, E)
            plt.plot(a.x1, a.v[4])
            plot.annotate_time(u_coarse.time)
            plt.savefig("{}/img/E-{:04d}.png".format(coarse_dir, k))

            plot = ScatterPlot()
            plot(coarse_grid, np.log10(rho))
            plot(coarse_grid, np.log10(rho0))
            plot.annotate_time(u_coarse.time)

            plt.plot(a.x1, np.log10(a.v[5]))
            plt.legend(["rho", "rho0", "rho_ref"])
            plt.savefig("{}/img/rho-{:04d}.png".format(coarse_dir, k))
            #
            # plot = ScatterPlot()
            # plot(coarse_grid, np.log10(p))
            # plot(coarse_grid, np.log10(p0))
            # plot.annotate_time(u_coarse.time)
            # plt.legend(["p", "p0"])
            # plt.savefig("{}/img/p-{:04d}.png".format(coarse_dir, k))

            plt.figure()
            radii = u_coarse.gravity["radii"]
            phi = u_coarse.gravity["phi"]
            phi = phi - phi[0]

            plt.plot(radii, phi, "+", label="phi")
            plt.plot(a.x1, a.v[23] - a.v[23][0], ".", label="phi_ref")

            text = "t = {:.3e}".format(u_coarse.time)
            plt.gca().text(0.70, 0.05, text, transform=plt.gca().transAxes)
            plt.legend()
            plt.savefig("{}/img/phi-{:04d}.png".format(coarse_dir, k))
            plt.close()

            plt.figure()
            radii = u_coarse.gravity["radii"][1:]
            phi = u_coarse.gravity["phi"][1:]
            phi = phi - phi[0]

            plt.loglog(radii, phi, "+", label="phi")
            # plt.loglog(a.x1, a.v[23] - a.v[23][0], ".", label="phi_ref")

            text = "t = {:.3e}".format(u_coarse.time)
            plt.gca().text(0.70, 0.05, text, transform=plt.gca().transAxes)
            plt.legend()
            plt.savefig("{}/img/phi-{:04d}.png".format(coarse_dir, k))
            plt.close()


class TableLabels:
    def __call__(self, col):
        return " ".join(str(v) for v in col.values())


def generate_grids():
    gmsh_template = read_txt("grids/janka.tmpl")

    for k, l in enumerate(mesh_levels):
        lc_min, lc_rel = mesh_sizes[k]
        geo = gmsh_template.replace("LC_MIN", str(lc_min))
        geo = geo.replace("RADIUS", str(polytrope_radius))
        geo = geo.replace("LC_REL", str(lc_rel))
        write_txt(geo_grid_name(l), geo)

    gmsh.generate_grids([geo_grid_name(l) for l in mesh_levels])


def main():
    parser = default_cli_parser("'janka' numerical experiment.")

    parser.add_argument(
        "--config-only",
        action="store_true",
        help="Write the config file to the working directory.",
    )

    args = parser.parse_args()

    if args.generate_grids:
        generate_grids()

    if args.run:
        build_zisa()
        # generate_grids()

        for c, r in all_runs:
            launch_all(c, force=args.force)

            if args.reference:
                launch_all(r, force=args.force)

    if args.restart:
        for c, r in all_runs:
            restart_all(c, args.restart_from)

            if args.reference:
                restart_all(r, args.restart_from)

    if args.post_process:
        for c, r in all_runs:
            post_process(c, r)

    if args.config_only:
        for c, r in all_runs:
            for cc in c:
                cc.save("config--" + cc.folder_name() + ".json")

            for rr in r:
                rr.save("config--" + rr.folder_name() + ".json")

    if args.copy_to_paper:
        dir = "${HOME}/git/papers/LucGrosheintz/papers/unstructured_well_balancing/img/janka"
        dir = os.path.expandvars(dir)

        for c, _ in all_runs:
            stem = c[0]["experiment"].short_id()
            files = glob.glob(stem + "*.tex")

            for f in files:
                shutil.copy(f, os.path.join(dir, os.path.basename(f)))


if __name__ == "__main__":
    main()
