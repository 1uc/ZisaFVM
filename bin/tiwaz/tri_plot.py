import matplotlib.pyplot as plt

def tri_plot(grid, data):
    print(grid.vertices.shape)
    print(grid.vertex_indices.shape)
    x, y = grid.vertices[:,0], grid.vertices[:,1]
    I = grid.vertex_indices

    ax = plt.gca()
    ax.tripcolor(x, y, I, facecolors=data)
