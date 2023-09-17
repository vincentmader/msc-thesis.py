from matplotlib.widgets import Slider
import numpy as np

from axis import DiscreteTimeAxis
from kernel import Kernel
from visualization.base import BasePlot


class EvolutionPlot(BasePlot):
    __slots__ = [
        "gs", "ax_1", "ax_2", "ax_3",
        "slider", "i_t", "t",
        "kernel", "N", "f", "m2f", "dm2f"
    ]

    def __init__(
        self, 
        kernel: Kernel,
        N: np.ndarray,
        f: np.ndarray,
        m2f: np.ndarray,
        dm2f: np.ndarray,
        *args, **kwargs
    ):
        super().__init__(
            figsize=(10, 8),
            *args, **kwargs
        )

        self.kernel = kernel
        self.N, self.f, self.m2f, self.dm2f = N, f, m2f, dm2f
        self.i_t = 0

        tg = DiscreteTimeAxis(kernel.cfg)
        tc = tg.grid_cell_centers
        self.t = tc

        self.gs = self.fig.add_gridspec(
            nrows=2, height_ratios=[2, 1], hspace=0.2,
            ncols=2, width_ratios=[20, 1], wspace=0.1,
        )
        self.ax_1 = self.fig.add_subplot(self.gs[0:1, 0:1])
        self.ax_2 = self.fig.add_subplot(self.gs[1:2, 0:1])
        self.ax_3 = self.fig.add_subplot(self.gs[0:2, 1:2])

        self.slider = Slider(
            ax=self.ax_3, orientation="vertical", label="t = 0", 
            valmin=0, valstep=1, valinit=0,
            valmax=self.N.shape[0] - 1, 
        )
        self.slider.on_changed(self.update)

    def draw(self):
        self.fig.canvas.draw_idle()

        mg = self.kernel.mg
        mb = mg.grid_cell_boundaries
        mc = mg.grid_cell_centers
        ac = mg.particle_radii

        n, N, M, dM = self.f, self.N, self.m2f, self.dm2f
        x, y = mc, M[self.i_t]

        self.ax_1.loglog(x, y, label="")
        self.ax_1.set_xlim(mb[0], mb[-1])
        self.ax_1.set_ylim(1e-12, 1e-8)
        self.ax_1.grid()
        self.ax_1.set_ylabel(r"dust particle density $\rho_i^s=m_i n_i \Delta m_i$")
        self.ax_1.set_title("temporal evolution of particle mass distribution")

        x, y = mc, dM[self.i_t]
        self.ax_2.loglog(x, y, label="")
        self.ax_2.set_xlim(mb[0], mb[-1])
        self.ax_2.set_ylim(1e-21, 1e-16)
        self.ax_2.grid()
        self.ax_2.set_xlabel("dust particle mass $m^c_i$")

    def update(self, i_t):
        self.i_t = i_t
        self.ax_1.clear()
        self.ax_2.clear()
        text = format_time(self.t[i_t])
        self.slider.label.set_text(f"t = {text}")
        self.draw()


def format_time(t) -> str:

    M = 60
    H = 60*M
    D = 24*H
    Y = 365.25*D
    KY = 1000 * Y
    MY = 1000 * KY

    if t >= MY:
        return f"{round(t / MY)} My"
    if t >= KY:
        return f"{round(t / KY)} ky"
    if t >= Y:
        return f"{round(t / Y)} y"
    if t >= D:
        return f"{round(t / D)} d"
    if t >= H:
        return f"{round(t / H)} h"
    if t >= M:
        return f"{round(t / M)} min"
    return f"{round(t)} s"
