#! /usr/bin/env python3

import matplotlib.pyplot as plt

from tiwaz.post_process import load_grid
from tiwaz.tri_plot import TriPlot

grid = load_grid(".")

plot = TriPlot()
plot.wire_frame(grid)
plt.gca().set_aspect("equal")

plot.save("grid.png")
