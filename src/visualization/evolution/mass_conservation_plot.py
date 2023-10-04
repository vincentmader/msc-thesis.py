import matplotlib.pyplot as plt

from constants import SECONDS_PER_YEAR
from visualization.base import BasePlot


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
        M = self.M
        err = (M - M[0]) / M[0]
        x, y = t, err

        # TODO Do with and without:
        # - Total mass error
        # - error per time
        # y = y[:-1] / (x[1:] - x[:-1])
        # x = x[:-1]

        plt.loglog(x, y, "red", label=r"$\Delta M_t>0$")
        plt.loglog(x, -y, "blue", label=r"$\Delta M_t<0$")
        plt.legend(loc="best")
        plt.xlabel("time $t$ [years]")
        plt.ylabel("mass error $\Delta M_t$ [kg]")
        plt.title(r"mass error $\Delta M_t=(M_t-M_0)/M_0$")

        color = "white" if self.cfg.mpl_dark_mode else "black" # TODO Make sure this works also when not calling from preset `p1`.
        rel_error = (M[-1]-M[0])/M[0]
        rel_error_sign = f"+" if rel_error > 0 else ""
        msgs = [
            ((0.5, 0.45), r"$M_i$ = " + f"{M[0]}"),
            ((0.5, 0.4), r"$M_f$ = " + f"{M[-1]}"),
            ((0.5, 0.35), r"$\frac{M_f-M_i}{M_i}$ = " + rel_error_sign + f"{rel_error*100} %"),
        ]
        for (x, y), msg in msgs:
            print(msg)
            plt_text(x, y, msg, color=color)


def plt_text(x, y, text, color="black"):
    plt.text(
        x, y, text,
        horizontalalignment="center",
        verticalalignment="center",
        transform=plt.gca().transAxes,
        color=color,
    )
