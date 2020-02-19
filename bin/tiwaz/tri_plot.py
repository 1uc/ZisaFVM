import numpy as np
import matplotlib.pyplot as plt


class TriPlot:
    def __init__(self):
        self.fig = plt.figure()

    def __del__(self):
        plt.close(self.fig)

    def wire_frame(self, grid):
        x, y = grid.vertices[:, 0], grid.vertices[:, 1]
        I = grid.vertex_indices

        plt.triplot(x, y, I, color="k", linewidth=0.1)

    def color_plot(self, grid, data):
        x, y = grid.vertices[:, 0], grid.vertices[:, 1]
        I = grid.vertex_indices

        masked_data = np.ma.masked_invalid(data)

        ax = plt.gca()
        tripcolor = ax.tripcolor(x, y, I, facecolors=masked_data)
        plt.colorbar(tripcolor)

    def quiver(self, grid, U, V):
        x, y = grid.cell_centers[:, 0], grid.cell_centers[:, 1]
        plt.quiver(x, y, U, V)

    def save(self, filename, dpi=300):
        plt.savefig(filename, dpi=dpi)

    def filename(self, data_filename, variable_key):
        return data_filename[:-3] + f"_tri-{variable_key}.png"


def tri_plot(grid, data):
    plot = TriPlot()
    plot.color_plot(grid, data)

    return plot
