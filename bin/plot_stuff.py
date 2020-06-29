import numpy as np
import matplotlib.pyplot as plt
from tiwaz.post_process import (
    load_grid,
    load_reference,
    load_data,
    find_last_data_file,
    find_steady_state_file,
)

cg = "gaussian_bump_3d_amp1.00e-01_o42222_RK4_isentropic_L0"
fg = "gaussian_bump_3d_amp1.00e-01_o42222_Fehlberg_isentropic_L3"

ref_grid = load_grid(f"{fg}/grids/gaussian_bump_3d-0/48/grid.h5")
approx_grid = load_grid(f"{cg}/grids/gaussian_bump_3d-0/2/grid.h5")

u_approx = load_data(find_last_data_file(cg), find_steady_state_file(cg))
u_ref = load_reference(fg, cg + "/gaussian_bump_3d-0")

plt.plot(np.linalg.norm(ref_grid.cell_centers, axis=1), u_ref.dvars["rho"], "k.")
plt.plot(np.linalg.norm(approx_grid.cell_centers, axis=1), u_approx.dvars["rho"], "r.")
plt.show()
