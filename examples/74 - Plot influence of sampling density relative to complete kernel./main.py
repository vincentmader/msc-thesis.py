from dataclasses import dataclass
import os, sys
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.widgets import Slider
import numpy as np
try:
    sys.path.append(os.path.join("..", "..", "src"))
    from models.axis import DiscreteMassAxis
    from models.axis import DiscreteTimeAxis
    from models.axis import AxisLabelVariant
    from config import Config
    from kernel import Kernel
    from utils.functions import format_time_as_years
    from visualization.base import BasePlot
    from visualization.preset import p1
except ModuleNotFoundError as e:
    raise e

Y_LIMITS = 1e-16, 1e5
COLORS   = [
    "red",
    "orange",
    "green",
    "blue",
    "violet",
    "black",
]
sampling_densities = [0.4, 0.6, 0.8, 1.0]

cfg = Config(
    enable_collision_sampling=True,
    # initial_mass_bin=40,
)
mg = DiscreteMassAxis(cfg)
tg = DiscreteTimeAxis(cfg)
mg = DiscreteMassAxis(cfg)
mb = mg.bin_boundaries
mc = mg.bin_centers
ac = mg.particle_radii

sampling_densities_and_numbers = [
    (rho_s, int(rho_s * cfg.mass_resolution**2 / 2))
    for rho_s in sampling_densities
]

m2fs, kernels, dm2fdts = {}, {}, {}
for rho_s, N_s in sampling_densities_and_numbers:
    if rho_s == 1.0:
        cfg.enable_collision_sampling = False
    else:
        cfg.enable_collision_sampling = True
        cfg.nr_of_samples = N_s
    kernel = Kernel(cfg)
    kernels[N_s] = kernel

    t, f, N, m2f, dm2fdt, M = p1.integrate(cfg, kernel)
    m2fs[N_s] = m2f
    dm2fdts[N_s] = dm2fdt

    print(f"sampling density {rho_s = }")
    print(f"max. nr. of sampled collisions {N_s = }")


@dataclass
class EvolutionMultiPlot(BasePlot):
    lines: dict[int, mpl.lines.Line2D]
    axes: dict[int, mpl.axis.Axis]
    i_t: int

    def __init__(self, *args, **kwargs):
        super().__init__(
            figsize=(16, 9),
            *args, **kwargs
        )

        self.axes, self.lines = {}, {}
        self.gs = self.fig.add_gridspec(
            nrows=1, height_ratios=[2], hspace=0.2,
            ncols=2, width_ratios=[30, 1], wspace=0.1,
        )
        self.ax_1 = self.fig.add_subplot(self.gs[0:1, 0:1])
        self.ax_3 = self.fig.add_subplot(self.gs[0:2, 1:2])
        self.ax_1c = self.ax_1.twiny()

        self.i_t = 100
        self.slider = Slider(
            ax=self.ax_3, orientation="vertical", label="t = 0", 
            valmin=0, valstep=1, valinit=self.i_t, valmax=tg.N-1, 
        )
        self.slider.on_changed(self.update)
        # self.ax_1.set_ylabel(r"accuracy error $\frac{\rho_{sampled} - \rho_{complete}}{\rho_{complete}}$")
        self.ax_1.set_ylabel(r"accuracy error $(\rho_{sampled} - \rho_{complete}) / \rho_{complete}$")
        # self.ax_1.set_xlabel("time $t$")
        self.ax_1.grid(True)

        plt.legend(loc="best")

    def draw(self):
        self.fig.canvas.draw_idle()

        for idx, (rho_s, N_s) in enumerate(sampling_densities_and_numbers):
            color = "black" if rho_s == 1.0 else COLORS[idx]

            cfg.nr_of_samples = N_s
            m2f = m2fs[N_s][self.i_t]
            dm2fdt = dm2fdts[N_s][self.i_t]

            m2f_complete = m2fs[int(mg.N**2 / 2)][self.i_t]
            y = (m2f - m2f_complete) / m2f_complete

            self.lines[2*N_s], = self.ax_1.loglog(
                mc, y, 
                label=r"$\rho_s=$" + f"{rho_s}",  # TODO
                color=color,
            )
            self.lines[2*N_s+1], = self.ax_1.loglog(
                mc, -y, "--", 
                color=color,
            )
            self.ax_1.set_xlim(mb[0], mb[-1])
            self.ax_1.set_ylim(Y_LIMITS[0], Y_LIMITS[1])
            self.ax_1.legend(loc="upper right")
            self.ax_1.grid(True)
            self.ax_1.set_xlabel("dust particle mass $m_i$ [m]")

            self.ax_1c.set_xscale("log")
            self.ax_1c.set_xlabel("dust particle radius $a_i$ [m]")
            self.ax_1c.xaxis.set_major_locator(ticker.LogLocator(base=10, numticks=15))
            self.ax_1c.set_xlim(ac[0], ac[-1])

    def update(self, i_t):
        self.i_t = i_t
        text = format_time_as_years(tg.bin_centers[i_t])
        self.slider.label.set_text(f"t = {text}")
        for idx, (rho_s, N_s) in enumerate(sampling_densities_and_numbers):

            m2f = m2fs[N_s][self.i_t]
            m2f_complete = m2fs[int(mg.N**2 / 2)][self.i_t]
            y = (m2f - m2f_complete) / m2f_complete

            self.lines[2*N_s].set_ydata(y)
            self.lines[2*N_s+1].set_ydata(-y)

p = EvolutionMultiPlot()
p.render()
