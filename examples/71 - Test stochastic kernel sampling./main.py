import os
import sys

import matplotlib.pyplot as plt
import numpy as np
try:
    sys.path.append(os.path.join("..", "..", "src"))
    from axis import DiscreteTimeAxis, DiscreteMassAxis
    from config import Config
    from disk import mass_distribution
    from disk.disk import disk_mass_from_distribution
    from kernel import Kernel
    from solver import Solver
    from visualization.v1.mass_error import DiskMassErrorPlot
    from visualization.v1.slider_plot_2 import InteractiveSliderLinePlot
except ModuleNotFoundError as e:
    raise e


def _run_integrator(kernel, K):
    mg = kernel.mg
    solver = Solver(cfg)
    n0 = mass_distribution.dirac_delta(cfg)
    N, f, m2f = solver.run(mg, n0, K)
    return N, f, m2f


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


if __name__ == "__main__":
    cfg = Config()
    kernel = Kernel(cfg)
    K, mc = kernel.K, kernel.mg.grid_cell_centers

    N, f, m2f = _run_integrator(kernel, K)

    mg = DiscreteMassAxis(cfg)
    dm = mg.grid_cell_widths
    tg = DiscreteTimeAxis(cfg)
    t = tg.grid_cell_centers
    Ms = [disk_mass_from_distribution(n, mc, dm) for n in f]

    # Calculate temporal derivative of mass distribution.
    dm2f = m2f[1:] - m2f[:-1]
    dm2f = list(dm2f)
    dm2f.append(dm2f[-1])  # TODO Fix array shapes in a better way than this.
    dm2f = np.array(dm2f)
    dm2f = [dm2f[i] / tg.grid_cell_widths[i] for i, _ in enumerate(dm2f)]

    plot_1(mc, m2f, dm2f)
    plot_2(t, Ms)
