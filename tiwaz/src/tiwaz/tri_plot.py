import numpy as np
import matplotlib.pyplot as plt


class TriPlot:
    def __init__(self, cmap="viridis", symmetric_limits=False):
        self.fig = plt.figure()
        self.cmap = cmap
        self.symmetric_limits = symmetric_limits

    def __del__(self):
        plt.close(self.fig)

    def wire_frame(self, grid):
        x, y = grid.vertices[:, 0], grid.vertices[:, 1]
        I = grid.vertex_indices

        plt.triplot(x, y, I, color="k", linewidth=0.1)

    def color_plot(self, grid, data, mask=True):
        x, y = grid.vertices[:, 0], grid.vertices[:, 1]
        I = grid.vertex_indices

        assert I.shape[1] == 3, "These can't be triangles."
        mask = np.logical_or(mask, np.logical_not(np.isreal(data)))
        masked_data = np.ma.masked_where(mask, data)

        if self.symmetric_limits:
            M = np.ma.max(np.ma.abs(masked_data))
            vmin, vmax = -M, M
        else:
            vmin, vmax = None, None

        ax = plt.gca()
        tripcolor = ax.tripcolor(
            x, y, I, facecolors=masked_data, cmap=self.cmap, vmin=vmin, vmax=vmax
        )

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
