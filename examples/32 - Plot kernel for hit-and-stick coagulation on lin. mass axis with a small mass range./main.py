import os
import sys

import numpy as np
try:
    sys.path.append(os.path.join("..", "..", "src"))
    from axis import DiscreteMassAxis
    from config import Config
    from dust import particle_radius_from_mass
    from kernel import Kernel
    from visualization.kernel.v3_2023_08_14.pcolor_matrix_subplot import PcolorMatrixSubplot
    from visualization.kernel.v3_2023_08_14.gridspec_plot import GridspecPlot
except ModuleNotFoundError as e:
    raise e


# Define kernel configuration.
cfg = Config(
    # Define mass axis.
    # On a linear grid, if we want to reach a mass grid spacing of exactly one,
    # we have to chose `mass_resolution = mass_max_value - mass_min_value`.
    mass_axis_scale="lin",
    mass_min_value=1,
    mass_max_value=13,
    mass_resolution=12,
    # Define processes to include in the simulation.
    enable_coagulation=True,
    enable_fragmentation=False,
    enable_cancellation_handling=False,
    enable_physical_collisions=False,
    relative_velocity_components=[],
)

# Define discrete mass axis.
mg = DiscreteMassAxis(cfg)

# Define kernel.
kernel = Kernel(cfg)
K = kernel.K

mg = kernel.mg
mc = mg.grid_cell_centers
rho_s = cfg.dust_particle_density
ac = particle_radius_from_mass(mc, rho_s)
N_m = mg.N

# Define list of kernels to plot.
kernels = [K]

i = np.arange(0, N_m, 1)

s1 = PcolorMatrixSubplot(
    i, i, kernel.K, 
    title="kernel gain contribution $G_{kij}$",
    xlabel="particle radius $a_j$ [m]",
    ylabel="particle radius $a_i$ [m]",
    scales=("lin", "lin", "lin"),
    symmetrized=True,
)

p = GridspecPlot([s1], add_slider=True)
p.render()
