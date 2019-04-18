import numpy as np
import matplotlib.pyplot as plt

class TriPlot:
    def __init__(self):
        self.fig = plt.figure()

    def color_plot(self, grid, data):
        x, y = grid.vertices[:,0], grid.vertices[:,1]
        I = grid.vertex_indices

        masked_data = np.ma.masked_invalid(data)
        print((min(masked_data), max(masked_data)))

        ax = plt.gca()
        tripcolor = ax.tripcolor(x, y, I, facecolors=masked_data, edgecolors="k")
        plt.colorbar(tripcolor)

    def quiver(self, grid, U, V):
        x, y = grid.cell_centers[:,0], grid.cell_centers[:,1]
        plt.quiver(x, y, U, V)


def tri_plot(grid, data):
    plot = TriPlot()
    plot.color_plot(grid, data)
