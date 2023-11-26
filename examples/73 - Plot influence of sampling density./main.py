from dataclasses import dataclass
import os, sys
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.widgets import Slider
import numpy as np
try:
    sys.path.append(os.path.join("..", "..", "src"))
    from axis import DiscreteMassAxis
    from axis import DiscreteTimeAxis
    from axis import AxisLabelVariant
    from config import Config
    from kernel import Kernel
    from visualization.base import BasePlot
    from visualization.preset import p1
except ModuleNotFoundError as e:
    raise e

# scale = mg.scale
# axis_label_variant = AxisLabelVariant.Radius if scale == "log" else AxisLabelVariant.Bin
# z_limits = (1e-20, 1e+10) if scale == "log" else (-1, 1)
# R = kernel.R_coag + kernel.R_frag

@dataclass
class EvolutionMultiPlot(BasePlot):
    lines: dict[int, mpl.lines.Line2D]
    axes: dict[int, mpl.axis.Axis]
    mg: DiscreteMassAxis
    tg: DiscreteTimeAxis
    i_t: int
    sampling_densities_and_numbers: list[tuple[float, int]]
    kernels: dict[int, Kernel]

    m2fs: dict[int, np.ndarray]

    def __init__(
        self, 
        # N: np.ndarray,
        # f: np.ndarray,
        # m2f: np.ndarray,
        # dm2fdt: np.ndarray,
        sampling_densities: list[float],
        *args, **kwargs
    ):
        super().__init__(
            figsize=(16, 9),
            *args, **kwargs
        )

        self.cfg = Config(
            enable_collision_sampling=True,
            # initial_mass_bin=40,
        )
        self.sampling_densities_and_numbers = [
            (rho_s, int(rho_s * self.cfg.mass_resolution**2))
            for rho_s in sampling_densities
        ]
        self.axes = {}
        self.lines = {}
        self.kernels = {}
        self.m2fs = {}

        self.mg = DiscreteMassAxis(self.cfg)
        self.tg = DiscreteTimeAxis(self.cfg)

        self.gs = self.fig.add_gridspec(
            nrows=2, height_ratios=[2, 1], hspace=0.2,
            ncols=2, width_ratios=[30, 1], wspace=0.1,
        )
        self.ax_1 = self.fig.add_subplot(self.gs[0:1, 0:1])
        self.ax_2 = self.fig.add_subplot(self.gs[1:2, 0:1])
        self.ax_3 = self.fig.add_subplot(self.gs[0:2, 1:2])
        self.ax_1c = self.ax_1.twiny()

        self.i_t = 190
        self.slider = Slider(
            ax=self.ax_3, orientation="vertical", label="t = 0", 
            valmin=0, valstep=1, valinit=self.i_t,
            valmax=self.tg.N-1, 
        )
        self.slider.on_changed(self.update)

        plt.grid(True)
        plt.ylabel(r"dust particle density $\rho_i^s=m_i n_i \Delta m_i$ [kg m$^{-3}$]")
        plt.legend(loc="best")

        for rho_s, N_s in self.sampling_densities_and_numbers:
            self.cfg.nr_of_samples = N_s
            if N_s == mg.N**2:
                pass
            kernel = Kernel(self.cfg)
            self.kernels[N_s] = kernel

            t, f, N, m2f, dm2fdt, M = p1.integrate(self.cfg, kernel)
            self.m2fs[N_s] = m2f

            print(f"sampling density {rho_s = }")
            print(f"max. nr. of sampled collisions {N_s = }")

    def draw(self):
        self.fig.canvas.draw_idle()

        mg = DiscreteMassAxis(self.cfg)
        mb = mg.bin_boundaries
        mc = mg.bin_centers
        ac = mg.particle_radii

        for rho_s, N_s in self.sampling_densities_and_numbers[::-1]:
            self.cfg.nr_of_samples = N_s
            m2f = self.m2fs[N_s][self.i_t]

            # print(self.m2fs[N_s])

            self.lines[2*N_s], = self.ax_1.loglog(
                mc, m2f, label=r"$\rho_s=$" + f"{rho_s}"
            )
            self.lines[2*N_s+1], = self.ax_1.loglog(mc, -m2f)
            # # self.lines_2, = self.ax_1.loglog(x, N[self.i_t], label=r"$n_i\Delta m_i=N_i$")
            # # self.lines_3, = self.ax_1.loglog(x, n[self.i_t], label=r"$n_i$") # TODO Plot these too?
            # self.ax_1.set_xlim(mb[0], mb[-1])
            # self.ax_1.set_ylim(1e-20, 1)  # TODO Generalize definition.
            self.ax_1.set_ylim(1e-16, 1e-9)
            self.ax_1.legend(loc="upper right")
            self.ax_1.grid(True)

            self.ax_1c.set_xscale("log")
            self.ax_1c.set_xlabel("dust particle radius $a_i$ [m]")
            # self.ax_1c.set_ylabel(r"$m_i \Delta n_i \Delta m_i$ [kg s$^{-2}$]") # TODO Why is this not shown?
            self.ax_1c.xaxis.set_major_locator(ticker.LogLocator(base=10, numticks=15))
            self.ax_1c.set_xlim(ac[0], ac[-1])

            # self.lines_4, = self.ax_2.loglog(mc, dM[self.i_t], label="> 0")
            # self.lines_5, = self.ax_2.loglog(mc, -dM[self.i_t], label="< 0")
            # self.ax_2.set_xlim(mb[0], mb[-1])
            # self.ax_2.set_ylim(1e-21, 1e-16)
            # self.ax_2.grid(True)
            # self.ax_2.set_xlabel("dust particle mass $m^c_i$ [kg]")
            # # self.ax_2.set_ylabel(r"$\frac{d\rho_i^s}{dt}$ [kg m$^{-3}$s$^{-1}$]")
            # self.ax_2.set_ylabel(r"d$\rho_i^s$ $/$ d$t$ [kg m$^{-3}$s$^{-1}$]")
            # self.ax_2.legend()

    def update(self, i_t):
        self.i_t = i_t
        # text = format_time(self.t[i_t])
        # self.slider.label.set_text(f"t = {text}")
        for rho_s, N_s in self.sampling_densities_and_numbers:
            self.lines[2*N_s].set_ydata(self.m2fs[N_s][i_t])
            self.lines[2*N_s+1].set_ydata(-self.m2fs[N_s][i_t])
        # # self.lines_5.set_ydata(self.m2f[i_t])
        # # self.lines_2.set_ydata(self.N[i_t])
        # # self.lines_3.set_ydata(self.f[i_t])
        # self.lines_4.set_ydata(self.dm2fdt[i_t])
        # self.lines_5.set_ydata(-self.dm2fdt[i_t])


sampling_densities = [0.3, 0.5, 0.7, 1.0]
# sampling_densities = [0.6, 1.0]
p = EvolutionMultiPlot(
    sampling_densities=sampling_densities,
)
p.render(close_plot=False, show_plot=False)
plt.show()
# a = lambda _: p.render(close_plot=False)
