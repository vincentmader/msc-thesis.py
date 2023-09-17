import os, sys
import numpy as np
try:
    sys.path.append(os.path.join("..", "..", "src"))
    from axis import DiscreteMassAxis, KernelAxis
    from config import Config
    from kernel import Kernel
    from visualization.base import GridspecPlot
    from visualization.kernel import KernelSubplot
except ModuleNotFoundError as e:
    raise e

# Define kernel configuration.
cfg = Config(
    # Define mass axis.
    # On a linear grid, if we want to reach a mass grid spacing of exactly one,
    # we have to chose `mass_resolution = mass_max_value - mass_min_value`.
    mass_axis_scale="lin",
    mass_min_value=1,
    mass_max_value=11,
    mass_resolution=10,
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
ac = mg.particle_radii
N_m = mg.N

# Define list of kernels to plot.
kernels = [K]

i = np.arange(0, N_m, 1)

s1 = KernelSubplot(
    mg, kernel.K,
    axis=KernelAxis.Bin,
    title="kernel $K_{kij}$",
    scales=("lin", "lin", "lin"),
    # symmetrized=True,
    cmap="bwr",
    z_limits=(-1, 1),
)

p = GridspecPlot([s1], add_slider=True)
p.render()
