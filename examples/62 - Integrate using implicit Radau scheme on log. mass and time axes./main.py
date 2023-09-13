import os
from pprint import pprint
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
    from visualization.v1.mass_error import DiskMassErrorPlot
    from visualization.v1.slider_plot_2 import InteractiveSliderLinePlot
except ModuleNotFoundError as e:
    raise e


# Load configuration from `../../config.toml`.
cfg = Config(enable_collision_sampling=False)

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

# ─────────────────────────────────────────────────────────────────────────────
# from disk.disk import Disk
# from disk.disk_region import DiskRegion
# from collision import collision_rate
# from kees_kernel import create_coag_kernel
# disk = Disk(cfg, rg, mg)
# disk_region = DiskRegion(cfg, disk)
# Cij = collision_rate(cfg, disk, disk_region)
# mgrain = mg.grid_cell_centers
# K = create_coag_kernel(mgrain, Cij)  # Kees
# ─────────────────────────────────────────────────────────────────────────────


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
        ylims_1=[1e-15, 1e-8],
        ylims_2=[1e-40, 1e-15],
    )
    plot.draw()
    plt.show()
    plt.close()


def plot_2(x, y):
    plot = DiskMassErrorPlot(x, y)
    plot.draw()
    plt.show()
    plt.close()


def plot_3(kernel, N, f, m2f, dm2f):
    from visualization.evolution.evolution import EvolutionPlot
    p = EvolutionPlot(kernel, N, f, m2f, dm2f)
    p.render()


if __name__ == "__main__":

    # Initialize mass distribution.
    n0 = mass_distribution.dirac_delta(cfg)
    # n0 = mass_distribution.mrn_distribution(cfg)

    # Run the solver.
    N, f, m2f = solver.run(mg, n0, K)

    # Calculate temporal derivative of mass distribution.
    dm2f = m2f[1:] - m2f[:-1]
    dm2f = list(dm2f)
    dm2f.append(dm2f[-1])  # TODO Fix array shapes in a better way than this.
    dm2f = np.array(dm2f)
    dm2f = [dm2f[i] / tg.grid_cell_widths[i] for i, _ in enumerate(dm2f)]

    # Prepare abscissa & ordinate for plot of disk mass error.
    t = tg.grid_cell_centers
    Ms = [disk_mass_from_distribution(n, mc, dm) for n in f]

    # Create plots.
    plot_3(kernel, N, f, m2f, dm2f)
    plot_1(mc, m2f, dm2f)
    plot_2(t, Ms)
