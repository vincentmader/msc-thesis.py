import os
import sys

import matplotlib.pyplot as plt
import numpy as np
try:
    sys.path.append(os.path.join("..", "..", "src"))
    from axis import DiscreteMassAxis, DiscreteRadialAxis, DiscreteTimeAxis
    from config import Config
    from disk import mass_distribution
    from disk.disk import disk_mass_from_distribution
    from kernel import Kernel
    from solver import Solver
    from visualization.evolution.v1.mass_error import DiskMassErrorPlot
    from visualization.evolution.v1.slider_plot_2 import InteractiveSliderLinePlot
except ModuleNotFoundError as e:
    raise e


# Define kernel configuration.
cfg = Config(
    mass_min_value=1,
    mass_max_value=51,
    mass_resolution=50,
    mass_axis_scale="lin",
    time_min_value=1,
    time_max_value=1e2,
    time_axis_scale="lin",
    solver_variant="explicit_euler",
    enable_coagulation=True,
    enable_fragmentation=False,
    enable_cancellation_handling=False,
    enable_physical_gas_density=False,
    enable_physical_collisions=False,
    relative_velocity_components=[],
    # ^ note: Enabling physical velocities leads to divergence.
)

# Define discrete axis for radial distance from star, as well as for mass.
rg = DiscreteRadialAxis(cfg)
mg = DiscreteMassAxis(cfg)
mc = mg.grid_cell_centers
dm = mg.grid_cell_widths

# Define kernel.
kernel = Kernel(cfg)
K = kernel.K

# Define temporal domain & solver.
tg = DiscreteTimeAxis(cfg)
solver = Solver(cfg)


def plot_1(m, m2f, dm2f):
    plot = InteractiveSliderLinePlot(
        cfg,
        m, m2f, dm2f,
        title="Temporal Evolution of Particle Mass Distribution Function",
        ylabel_1="$n_i\cdot m_i\cdot\Delta m_i$    [kg m $^{-3}$]",
        xlabel_1="mass $m_i$ [kg]",
        ylabel_2=r"$\frac{\Delta n_i\cdot m_i\cdot\Delta m_i}{\Delta t}$   [kg m$^{-3}$ s$^{-1}$]",
        xlabel_2="mass $m_i$ [kg]",
        xlims_1=(m[0], m[-1]),
        xlims_2=(m[0], m[-1]),
    )
    plot.draw()
    plt.show()
    plt.close()


def plot_2(x, y):
    plot = DiskMassErrorPlot(x, y)
    plot.draw()
    plt.show()
    plt.close()


if __name__ == "__main__":

    # Initialize mass distribution.
    n0 = mass_distribution.dirac_delta(cfg)

    # Run the solver.
    N, f, m2f = solver.run(mg, n0, K)

    # Calculate temporal derivative of mass distribution.
    dm2f = m2f[1:] - m2f[:-1]
    dm2f = list(dm2f)
    dm2f.append(dm2f[-1])  # TODO Fix array shapes in a better way than this.
    dm2f = np.array(dm2f)

    # Prepare abscissa & ordinate for plot of disk mass error.
    t = tg.grid_cell_centers
    Ms = [disk_mass_from_distribution(n, mc, dm) for n in f]

    # Create plots.
    plot_1(mc, m2f, dm2f)
    plot_2(t, Ms)