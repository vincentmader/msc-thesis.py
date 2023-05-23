import os
import sys

import numpy as np
try:
    sys.path.append(os.path.join("..", "..", "src"))
    from config import Config
    from disk import MassGrid, Disk, RadialGrid, mass_distribution, DiskRegion, TimeGrid
    from disk.disk import disk_mass_from_distribution
    from kernel import Kernel
    from kernel import collision_rate
    from kernel.collision_rate import collision_rate
    from solver import Solver
    from utils.plotting import plt_show_then_close
    from visualization.mass_error import DiskMassErrorPlot
    from visualization.slider_plot_2 import InteractiveSliderLinePlot
except ModuleNotFoundError as e:
    raise e


# Define kernel configuration.
cfg = Config(
    mass_min_value=1,
    mass_max_value=1e6,
    mass_resolution=50,
    mass_axis_scale="log",
    time_min_value=1,
    time_max_value=1e2,
    time_axis_scale="lin",
    solver_variant="explicit_euler",
    enable_coagulation=True,
    enable_fragmentation=False,
    enable_cancellation_handling=False,
    enable_physical_gas_density=False,
    enable_physical_cross_sections=False,
    enable_physical_relative_velocities=[],
    # ^ NOTE: Enabling physical velocities leads to divergence.
)

# Define discrete axis for radial distance from star, as well as for mass.
rg = RadialGrid(cfg)
mg = MassGrid(cfg)
m = mg.grid_cell_boundaries()[:-1]  # TODO use bounds or centers?

# Define disk, the position of interest in it, & the disk properties there.
disk = Disk(cfg, rg, mg)
disk_region = DiskRegion(cfg, disk)
R_coll = collision_rate(cfg, disk, disk_region) # TODO Define dynamically based on cfg.

# Define kernel.
kernel = Kernel(cfg)
K = kernel.K(mg, R_coll)

# Define temporal domain & solver.
tg = TimeGrid(cfg)
solver = Solver(cfg)


def plot_1(m, m2f, dm2f):
    plot = InteractiveSliderLinePlot(
        cfg,
        m, m2f, dm2f,
        title="Temporal Evolution of Particle Mass Distribution Function",
        ylabel_1="$m_i^2 n_i$ [kg m $^{-3}$]",
        xlabel_1="mass $m_i$ [kg]",
        ylabel_2=r"$\frac{dn_i}{dt}\cdot\Delta t$",
        xlabel_2="mass $m_i$ [kg]",
        xlims_1=(m[0], m[-1]),
        xlims_2=(m[0], m[-1]),
    )
    plot.draw()
    plt_show_then_close()


def plot_2(x, y):
    plot = DiskMassErrorPlot(x, y)
    plot.draw()
    plt_show_then_close()


if __name__ == "__main__":

    # Initialize mass distribution.
    n0 = mass_distribution.dirac_delta(cfg)

    # Run the solver.
    N, f, m2f = solver.run(mg, n0, K)

    # Calculate temporal derivative of mass distribution.
    dm2f = m2f[1:] - m2f[:-1]
    dm2f = list(dm2f) 
    dm2f.append(dm2f[-1]) # TODO Fix array shapes in a better way than this.
    dm2f = np.array(dm2f)

    # Prepare abscissa & ordinate for plot of disk mass error.
    t = tg.grid_cell_centers()
    Ms = [disk_mass_from_distribution(m, n) for n in f]

    # Create plots.
    plot_1(m, m2f, dm2f)
    plot_2(t, Ms)
