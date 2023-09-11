import os
import sys
import numpy as np
try:
    sys.path.append(os.path.join("..", "..", "src"))
    from config import Config
    from kernel import Kernel
    from visualization.kernel.v3_2023_08_14.pcolor_matrix_subplot import PcolorMatrixSubplot
    from visualization.kernel.v3_2023_08_14.gridspec_plot import GridspecPlot
except ModuleNotFoundError as e:
    raise e


cfg = Config(
    mass_axis_scale="lin",
    mass_min_value=1,
    mass_max_value=50,
    enable_fragmentation=False,
    enable_physical_collisions=False,
    relative_velocity_components=[],
)
kernel_1 = Kernel(cfg)

ac = np.linspace(0, 49, 50)

s1 = PcolorMatrixSubplot(
    ac, ac, kernel_1.K_gain, 
    title="kernel gain contribution $G_{kij}$",
    xlabel="particle radius $a_j$ [m]",
    ylabel="particle radius $a_i$ [m]",
    scales=("lin", "lin", "lin"),
    symmetrize=True,
)
s2 = PcolorMatrixSubplot(
    ac, ac, -kernel_1.K_loss,
    title="kernel loss contribution $L_{kij}$",
    xlabel="particle radius $a_j$ [m]",
    scales=("lin", "lin", "lin"),
    symmetrize=True,
)
subplots = [s1, s2]

p = GridspecPlot(subplots, add_slider=True)
p.render()
