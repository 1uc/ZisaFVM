import os
import glob

import numpy as np
import h5py

from .launch_params import folder_name
from .xdmf import generate_xdmf, xml_to_string
from .io import read_json


class Grid:
    def __init__(self, grid_name):
        print(grid_name)
        with h5py.File(grid_name, "r") as h5:
            self.vertex_indices = np.array(h5["vertex_indices"])
            self.vertices = np.array(h5["vertices"])

            if "cell_flags/ghost_cell" in h5:
                self.is_ghost_cell = np.array(h5["cell_flags/ghost_cell"]).astype(bool)

            if "cell_centers" in h5:
                self.cell_centers = np.array(h5["cell_centers"])

            if "volumes" in h5:
                self.volumes = np.array(h5["volumes"])
                self.volumes = np.ma.masked_array(self.volumes, mask=self.is_ghost_cell)

            if "dx_max" in h5:
                self.dx_max = h5["dx_max"][()]

            self.n_cells = self.vertex_indices.shape[0]


class Snapshot:
    def __init__(self, data_filename, delta_filename=None, steady_state_filename=None):

        cvar_keys = ["rho", "mv1", "mv2", "E"]
        avar_keys = [f"mq{i}" for i in range(128)]
        xvar_keys = [
            "v1",
            "v2",
            "v3",
            "p",
            "cs",
            "h",
            "s",
            "E_th",
            "E_p",
            "p_p",
            "p_th",
        ]

        max_avars = 8
        xvar_keys = (
            xvar_keys
            + [f"mq{i}" for i in range(max_avars)]
            + [f"q{i}" for i in range(max_avars)]
        )

        self.cvars = dict()
        self.xvars = dict()
        self.dvars = dict()

        with h5py.File(data_filename, "r") as h5:
            if "time" in h5:
                self.time = h5["time"][()]

            for key in cvar_keys:
                self.cvars[key] = np.array(h5[key])

            for key in xvar_keys:
                if key in h5:
                    self.xvars[key] = np.array(h5[key])

            self.gravity = dict()
            if "model/gravity/radii" in h5:
                self.gravity["radii"] = np.array(h5["model/gravity/radii"])

            if "model/gravity/phi" in h5:
                self.gravity["phi"] = np.array(h5["model/gravity/phi"])

        if delta_filename:
            with h5py.File(delta_filename, "r") as h5:
                for key in cvar_keys:
                    self.dvars[key] = np.array(h5[key])

        elif steady_state_filename:
            with h5py.File(steady_state_filename, "r") as h5:
                for key in cvar_keys:
                    self.dvars[key] = self.cvars[key] - np.array(h5[key])

    def __getitem__(self, key):
        if key in self.cvars:
            return self.cvars[key]
        elif key in self.xvars:
            return self.xvars[key]
        elif key[1:] in self.dvars:
            return self.dvars[key[1:]]
        else:
            raise Exception(f"Unknown variable. [{key}]")


def load_grid(path):
    if os.path.isfile(path):
        return Grid(path)
    elif os.path.isdir(path):
        return Grid(find_grid(path))
    else:
        raise Exception(f"Invalid path. [{path}]")


def load_reference(reference_dir, grid_name):
    stem, _ = os.path.splitext(os.path.basename(grid_name))
    down_dir = os.path.join(reference_dir, "down_sampled", stem)
    ref_filename = os.path.join(down_dir, "reference.h5")
    delta_filename = os.path.join(down_dir, "delta.h5")
    return Snapshot(ref_filename, delta_filename=delta_filename)


def load_data(data_name, steady_state_filename=None):
    return Snapshot(data_name, steady_state_filename=steady_state_filename)


def find_grid(directory):
    f = os.path.join(directory, "grid.h5")
    if os.path.isfile(f):
        return f

    c = os.path.join(directory, "config.json")
    if os.path.isfile(c):
        config = read_json(c)
        g = os.path.join(directory, config["grid"]["file"])

        if os.path.isfile(g):
            return g

        g = next(iter(glob.glob(os.path.join(g, "grid.h5"))))
        if os.path.isfile(g):
            return g

    raise Exception("Failed to find the grid.")


def find_data_files(directory):
    return sorted(glob.glob(os.path.join(directory, "*_data-*.h5")))


def find_last_data_file(directory):
    return max(find_data_files(directory))


def find_first_data_file(directory):
    return min(find_data_files(directory))


def find_steady_state_file(directory):
    return os.path.join(directory, "steady_state.h5")


def extract_solver_data(solver, data):
    """Returns the sub sequence of data for one solver.

    This is useful if one needs to plot the data for one solver
    across all grid levels.

    Assume the data is a list of dictionaries, e.g.
      [
        {"order": 1, "wb": "isentropic", "grid_level": 2, "bla": 42.0, ...},
        {"order": 1, "wb": "isentropic", "grid_level": 1, "bla": 42.0, ...},
        {"order": 2, "wb": "isentropic", "grid_level": 1, "bla": 42.0, ...},
        {"order": 2, "wb": "isentropic", "grid_level": 2, "bla": 42.0, ...}
      ]

    Then `extract_solver_data({"order": 1, "wb": "isentropic"})` will return
      [
        {"order": 1, "wb": "isentropic", "grid_level": 1, "bla": 42.0, ...},
        {"order": 1, "wb": "isentropic", "grid_level": 2, "bla": 42.0, ...}
      ]

    """
    filtered_data = filter(lambda r: all(r[key] == solver[key] for key in solver), data)
    return sorted(filtered_data, key=lambda r: r["grid_level"])


def load_results(coarse_runs, reference_run):
    results = list()

    ref_dir = folder_name(reference_run)
    fine_grid = load_grid(ref_dir)

    u_exact = load_data(find_last_data_file(ref_dir), find_steady_state_file(ref_dir))
    t_end = u_exact.time

    for coarse_run in coarse_runs:
        coarse_dir = folder_name(coarse_run)

        coarse_grid = load_grid(coarse_dir)
        volumes = coarse_grid.volumes

        u_approx = load_data(
            find_last_data_file(coarse_dir), find_steady_state_file(coarse_dir)
        )
        u_ref = load_reference(ref_dir, coarse_run.grid_filename())

        rho_approx = u_approx.cvars["rho"]
        drho_approx = u_approx.dvars["rho"]

        rho_ref = u_ref.cvars["rho"]
        drho_ref = u_ref.dvars["rho"]

        l1_err = np.sum(volumes * np.abs(rho_approx - rho_ref))
        l1_eq_err = np.sum(volumes * np.abs(drho_approx - drho_ref))

        # if np.abs(u_approx.time - t_end) > 1e-8 * t_end:
        #     l1_err = np.nan * l1_err
        #     l1_eq_err = np.nan * l1_eq_err

        results.append(
            {
                "order": coarse_run.order(),
                "wb": coarse_run.well_balancing(),
                "grid_level": coarse_run.grid_level(),
                "dx_max": coarse_grid.dx_max,
                "grid": coarse_grid,
                "fine_grid": fine_grid,
                "u_approx": u_approx,
                "u_ref": u_ref,
                "u_exact": u_exact,
                "l1_error": l1_err,
                "l1_eq_error": l1_eq_err,
            }
        )

        print(results[-1])

    orders = sorted(set(r["order"] for r in results))
    wbs = sorted(set(r["wb"] for r in results))

    columns = []
    for o in orders:
        for wb in wbs:
            columns.append({"order": o, "wb": wb})

    return results, columns


def find_components(data_file):
    keys = ["rho", "mv1", "mv2", "mv3", "E", "cs", "h", "p", "s"]
    with h5py.File(data_file, "r") as h5:
        return [k for k in keys if k in h5]


def write_xdmf(scheme):
    dirname = folder_name(scheme)
    grid_name = find_grid(dirname)
    xdmf_name = f"{dirname}/{scheme['experiment']['name']}.xdmf"

    data_files = find_data_files(dirname)
    components = find_components(data_files[0])

    xdmf_str = xml_to_string(generate_xdmf(grid_name, data_files, components))
    with open(xdmf_name, "w") as f:
        f.write(xdmf_str)
