import os
import sys

import numpy as np
try:
    sys.path.append(os.path.join("..", "..", "src"))
    from config import Config
    from dust import particle_radius_from_mass
    from kernel import Kernel
    from visualization.kernel.v3_2023_08_14.pcolor_matrix_subplot import PcolorMatrixSubplot
    from visualization.kernel.v3_2023_08_14.gridspec_plot import GridspecPlot
except ModuleNotFoundError as e:
    raise e

cfg_1 = Config(
    enable_cancellation_handling=False,
)
kernel_1 = Kernel(cfg_1)
K_1 = kernel_1.K

cfg_2 = Config(
    enable_cancellation_handling=True,
)
kernel_2 = Kernel(cfg_2)
K_2 = kernel_2.K

mg = kernel_1.mg
mc = mg.grid_cell_centers
rho_s = cfg_1.dust_particle_density
ac = particle_radius_from_mass(mc, rho_s)

K_diff = K_1 - K_2
K_equal = K_diff < 1e-14

s1 = PcolorMatrixSubplot(
    ac, ac, K_1,
    title="cancellation handling deactivated",
    xlabel="particle radius $a_j$ [m]",
    ylabel="particle radius $a_i$ [m]",
    scales=("log", "log", "lin"),
    symmetrized=True,
)
s2 = PcolorMatrixSubplot(
    ac, ac, K_2,
    title="cancellation handling activated",
    xlabel="particle radius $a_j$ [m]",
    scales=("log", "log", "lin"),
    symmetrized=True,
)
s3 = PcolorMatrixSubplot(
    ac, ac, np.abs(K_diff),
    title="abs($K_{canc}-K_{nocanc}$)",
    xlabel="particle radius $a_j$ [m]",
    scales=("log", "log", "log"),
    symmetrized=True,
)

p = GridspecPlot([s1, s2, s3], add_slider=True)
p.render()

# ═════════════════════════════════════════════════════════════════════════════

from visualization.kernel.v2.mass_conservation import KernelMassConservationPlot

# Plot `\sum_{ij} m_k \Delta m_k K_kij`.
p = KernelMassConservationPlot(cfg_1, mg, K_1)
p.show()

# Plot `\sum_{ij} m_k \Delta m_k K_kij`.
p = KernelMassConservationPlot(cfg_2, mg, K_2)
p.show()
