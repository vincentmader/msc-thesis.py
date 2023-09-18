import matplotlib.pyplot as plt

from visualization.base import BasePlot


class MassConservationPlot(BasePlot):
    __slots__ = [
        "gs", "ax_1", "ax_2", 
        "t", "M",
    ]

    def __init__(self, t, M, *args, **kwargs):
        self.t = t
        self.M = M
        super().__init__(*args, **kwargs)

    def draw(self):
        t = self.t / (365 * 24 * 60 * 60)
        M = self.M

        err = (M - M[0]) / M[0]

        x = t
        y = err

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

        rel_error = (M[-1]-M[0])/M[0]
        rel_error_sign = f"+" if rel_error > 0 else ""
        msgs = [
            ((0.5, 0.45), r"$M_i$ = " + f"{M[0]}"),
            ((0.5, 0.4), r"$M_f$ = " + f"{M[-1]}"),
            ((0.5, 0.35), r"$\frac{M_f-M_i}{M_i}$ = " + rel_error_sign + f"{rel_error*100} %"),
        ]
        for (x, y), msg in msgs:
            print(msg)
            plt_text(x, y, msg)


def plt_text(x, y, text):
    plt.text(
        x, y, text,
        horizontalalignment='center',
        verticalalignment='center',
        transform=plt.gca().transAxes,
    )