import numpy as np
import matplotlib.pyplot as plt

from .post_process import extract_solver_data
from .colors import graded_colors


class ScatterPlot:
    def __init__(self):
        self.fig = plt.figure()

    def __del__(self):
        plt.close(fig=self.fig)

    def __call__(self, grid, data, color=None):
        radii = np.linalg.norm(grid.cell_centers, axis=1)

        plt.figure(self.fig.number)
        plt.plot(radii, data, "*", markersize=2, color=color)

    def reference(self, grid, data):
        radii = np.linalg.norm(grid.cell_centers, axis=1)

        plt.figure(self.fig.number)
        plt.plot(radii[::100], data[::100], "k.", markersize=2)

    def save(self, filename):
        plt.figure(self.fig.number)
        plt.savefig(filename)

    def filename(self, data_filename, variable_key):
        return data_filename[:-3] + f"_scatter-{variable_key}.png"

    def finalize(self, label):
        plt.figure(self.fig.number)
        plt.title(label)

    def annotate_time(self, time):
        plt.figure(self.fig.number)
        ax = plt.gca()

        text = "t = {:.3e}".format(time)
        ax.text(0.05, 0.05, text, transform=ax.transAxes)


def scatter_plot(grid, data):
    plot = ScatterPlot()
    plot(grid, data)

    return plot


def plot_visual_convergence(data, solvers, labels, filename):
    for solver in solvers:
        sdata = extract_solver_data(solver, data)
        colors = graded_colors("blue", len(sdata))

        for var in ["rho", "E"]:
            plot = ScatterPlot()

            for d, c in zip(sdata, colors):
                plot(d["grid"], d["u_approx"].dvars[var], color=c)
            plot.reference(d["fine_grid"], d["u_exact"].dvars[var])

            label = labels(solver)
            plot.finalize(label)
            plot.save(filename + "_" + label.replace(" ", "_") + "_" + var + ".png")
