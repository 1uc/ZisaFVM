import numpy as np
import matplotlib.pyplot as plt

def tri_plot(grid, data):
    print(grid.vertices.shape)
    print(grid.vertex_indices.shape)
    x, y = grid.vertices[:,0], grid.vertices[:,1]
    I = grid.vertex_indices

    masked_data = np.ma.masked_invalid(data)
    print((min(masked_data), max(masked_data)))

    ax = plt.gca()
    tripcolor = ax.tripcolor(x, y, I, facecolors=masked_data)
    plt.colorbar(tripcolor)
