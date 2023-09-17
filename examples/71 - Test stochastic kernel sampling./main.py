import os, sys
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
    from visualization.evolution import EvolutionPlot, MassConservationPlot
    from visualization.v1.mass_error import DiskMassErrorPlot
    from visualization.v1.slider_plot_2 import InteractiveSliderLinePlot
except ModuleNotFoundError as e:
    raise e

def _run_integrator(kernel, K):
    mg = kernel.mg
    solver = Solver(cfg)
    n0 = mass_distribution.dirac_delta(cfg)
    N, f, m2f, dm2f = solver.run(n0, K)
    return N, f, m2f, dm2f


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
    p = EvolutionPlot(kernel, N, f, m2f, dm2f)
    p.render()


def plot_4(t, Ms):
    p = MassConservationPlot(t, Ms)
    p.render()


if __name__ == "__main__":
    cfg = Config()
    kernel = Kernel(cfg)
    K, mc = kernel.K, kernel.mg.bin_centers

    N, f, m2f, dm2f = _run_integrator(kernel, K)

    mg = DiscreteMassAxis(cfg)
    dm = mg.bin_widths
    tg = DiscreteTimeAxis(cfg)
    t = tg.bin_centers
    Ms = [disk_mass_from_distribution(n, mc, dm) for n in f]

    # plot_1(mc, m2f, dm2f)
    # plot_2(t, Ms)
    plot_3(kernel, N, f, m2f, dm2f)
    plot_4(t, Ms)
