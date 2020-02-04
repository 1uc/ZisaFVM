import numpy as np
import matplotlib.pyplot as plt

from .post_process import extract_solver_data


class ScatterPlot:
    def __init__(self):
        self.fig = plt.figure()

    def __del__(self):
        plt.close(fig=self.fig)

    def __call__(self, grid, data):
        radii = np.linalg.norm(grid.cell_centers, axis=1)

        plt.figure(self.fig.number)
        plt.plot(radii, data, "r*", markersize=2)

    def reference(self, grid, data):
        radii = np.linalg.norm(grid.cell_centers, axis=1)

        plt.figure(self.fig.number)
        plt.plot(radii[::100], data[::100], "k.", markersize=2)

    def save(self, filename):
        plt.figure(self.fig.number)
        plt.savefig(filename)

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

        for var in ["rho", "E"]:
            plot = ScatterPlot()

            for d in sdata:
                plot(d["grid"], d["u_approx"].dvars[var])
            # plot.reference(sdata[0]["fine_grid"], sdata[0]["u_exact"].dvars[var])

            label = labels(solver)
            plot.finalize(label)
            plot.save(filename + "_" + label.replace(" ", "_") + "_" + var + ".png")
