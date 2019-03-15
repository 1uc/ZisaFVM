import numpy as np
import matplotlib.pyplot as plt

from . post_process import extract_solver_data

class ScatterPlot:
    def __init__(self):
        self.fig = plt.figure()

    def __del__(self):
        plt.close(fig=self.fig)

    def __call__(self, grid, data):
        radii = np.linalg.norm(grid.cell_centers, axis=1)
        plt.plot(radii, data, "*", markersize=2)

    def reference(self, grid, data):
        radii = np.linalg.norm(grid.cell_centers, axis=1)

        plt.figure(self.fig.number)
        plt.plot(radii[::100], data[::100], 'k.', markersize=2)

    def save(self, filename):
        plt.figure(self.fig.number)
        plt.savefig(filename)

    def finalize(self, label):
        plt.figure(self.fig.number)
        plt.title(label)


def plot_visual_convergence(data, solvers, labels, filename):
    for solver in solvers:
        sdata = extract_solver_data(solver, data)

        plot = ScatterPlot()

        for d in sdata:
            plot(d["grid"], d["u_approx"].dvars["rho"])

        plot.reference(d["fine_grid"], d["u_exact"].dvars["rho"])

        label = labels(solver)
        plot.finalize(label)
        plot.save(filename + "_" + label + ".png")
