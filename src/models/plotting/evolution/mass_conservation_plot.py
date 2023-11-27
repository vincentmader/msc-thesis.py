import matplotlib.pyplot as plt

from constants import SECONDS_PER_YEAR
from models.plotting.base import BasePlot


class MassConservationPlot(BasePlot):
    __slots__ = [
        "gs", "ax_1", "ax_2", 
        "cfg", "t", "M",
    ]

    def __init__(self, cfg, t, M, *args, **kwargs):
        self.cfg = cfg
        self.t = t
        self.M = M
        super().__init__(*args, **kwargs)

    def draw(self):
        t = self.t / SECONDS_PER_YEAR
        M = self.M  # TODO Rename: M -> rho
        err = (M - M[0]) / M[0] * 100
        x, y = t, err

        # TODO Do with and without:
        # - Total mass error
        # - error per time
        # y = y[:-1] / (x[1:] - x[:-1])
        # x = x[:-1]

        plt.loglog(x,  y, "red",  label=r"$\Delta \rho_t>0$")
        plt.loglog(x, -y, "blue", label=r"$\Delta \rho_t<0$")
        plt.legend(loc="best")
        plt.xlabel("time $t$ [years]")
        plt.ylabel(r"mass error $\Delta\rho_t$ [%]")
        plt.grid()

        rel_error = (M[-1] - M[0]) / M[0]
        rel_error_sign = f"+" if rel_error > 0 else ""
        title = r"mass error $\Delta\rho_t=(\rho_t-\rho_0)/\rho_0$ = " + rel_error_sign + f"{rel_error*100} %"
        plt.title(title)


def plt_text(x, y, text, color="black"):
    plt.text(
        x, y, text,
        horizontalalignment="center",
        verticalalignment="center",
        transform=plt.gca().transAxes,
        color=color,
    )
