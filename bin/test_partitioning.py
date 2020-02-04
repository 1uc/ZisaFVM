#! /usr/bin/env python3
import numpy as np
import h5py

import matplotlib.pyplot as plt

import tiwaz.xdmf as xdmf
from tiwaz.tri_plot import TriPlot
from tiwaz.post_process import load_grid

# -- test the partitioned grid.
grid_file = "grid_full.h5"
data_files = ["f0.h5"]

with h5py.File(grid_file, "a") as h5:
    vertex_indices = np.array(h5["vertex_indices"])
    vertices = np.array(h5["vertices"])

    n_cells = vertex_indices.shape[0]
    n_vertices = vertices.shape[0]

    if not "n_cells" in h5:
        h5["n_cells"] = n_cells
    if not "n_vertices" in h5:
        h5["n_vertices"] = n_vertices
    if not "max_neighbours" in h5:
        h5["max_neighbours"] = vertex_indices.shape[1]


components = ["rho", "E"]
with open("test.xdmf", "w") as f:
    f.write(xdmf.xml_to_string(xdmf.generate_xdmf(grid_file, data_files, components)))

# test the subgrid
grid_files = [f"grid_part-{n:04d}.h5" for n in range(12)]
data_files = [f"data_part-{n:04d}.h5" for n in range(12)]

for grid_file, data_file in zip(grid_files, data_files):
    grid = load_grid(grid_file)
    with h5py.File(data_file, "r") as h5:
        data = np.array(h5["upsilon"])
        upsilon = data[:, 0]

    plot = TriPlot()
    plot.wire_frame(grid)
    plot.color_plot(grid, upsilon)
    plt.show()
