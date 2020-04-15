#! /usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


import h5py
import argparse

import tiwaz
from tiwaz.post_process import find_grid
from tiwaz.post_process import find_steady_state_file
from tiwaz.post_process import find_data_files
from tiwaz.post_process import load_grid


def load(data_file):
    steady_state_file = find_steady_state_file(os.path.dirname(data_file))
    if not os.path.isfile(steady_state_file):
        steady_state_file = None

    return tiwaz.post_process.load_data(data_file, steady_state_file)


def load_array(h5_file, key):
    with h5py.File(h5_file, "r") as h5:
        return np.array(h5[key])

def draw_line(xi, xj, **plot_kwargs):
    ax = plt.gca()
    xij = np.array([xi, xj])
    ax.plot(xij[:,0], xij[:,1], xij[:,2], **plot_kwargs)

def plot_stencil(grid, indices, color="red"):
    i = indices[0]
    xi = grid.cell_centers[i]

    ax = plt.gca()

    for j in indices[1:]:
        xj = grid.cell_centers[j]
        draw_line(xi, xj, color=color)


def plot_wireframe(grid, indices):
    ax = plt.gca()
    vi = grid.vertex_indices
    v = grid.vertices

    for i in indices:
        for k in vi[i,:]:
            for l in vi[i,:]:
                if k != l:
                    draw_line(v[k], v[l], color="black")


def write_stencil_mask(outname, grid, stencils):
    n_cells = grid.n_cells
    stencil_masks = np.zeros((n_cells, len(stencils)))

    for k, s in enumerate(stencils):
        stencil_masks[s,k] = 1.0
        stencil_masks[s[0], k] = 2.0

    with h5py.File(outname, "w") as h5:
        for k in range(stencil_masks.shape[1]):
            h5[f"s{k}"] = stencil_masks[:,k]


if __name__ == "__main__":
    parser_help = "Visualize the output Zisa."
    parser = argparse.ArgumentParser(description=parser_help)

    parser.add_argument(
        "--grid", required=True, help="The filename of the grid.",
    )

    parser.add_argument(
        "--index", type=int, required=True, help="Index of the cell whose stencils you want.",
    )



    args = parser.parse_args()

    grid_name = args.grid
    grid = load_grid(grid_name)


    # fig = plt.figure()
    # ax = fig.add_subplot(111, projection='3d')

    # indices = load_array(grid_name, "stencils/0/0")
    # plot_stencil(grid, indices)
    # plot_wireframe(grid, indices)

    # indices = load_array(grid_name, "stencils/0/1")
    # plot_stencil(grid, indices, color="green")
    # plot_wireframe(grid, indices)

    # indices = load_array(grid_name, "stencils/0/2")
    # plot_stencil(grid, indices, color="blue")
    # plot_wireframe(grid, indices)

    # indices = load_array(grid_name, "stencils/0/3")
    # plot_stencil(grid, indices, color="orange")
    # plot_wireframe(grid, indices)

    # indices = load_array(grid_name, "stencils/0/4")
    # plot_stencil(grid, indices, color="purple")
    # plot_wireframe(grid, indices)

    # plt.show()

    stencils = [load_array(grid_name, "stencils/{}/{}".format(args.index, k)) for k in range(5)]
    write_stencil_mask("__unit_tests--stencil_data.h5", grid, stencils)
