import matplotlib.pyplot as plt

from visualization.base import BasePlot


class MassConservationPlot(BasePlot):

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
        y = y[:-1] / (x[1:] - x[:-1])
        x = x[:-1]

        plt.loglog(x, y, "red", label=r"$\Delta M_t>0$")
        plt.loglog(x, -y, "blue", label=r"$\Delta M_t<0$")
        plt.legend(loc="best")
        plt.xlabel("time $t$ [years]")
        plt.ylabel("mass error $\Delta M_t$ [kg]")

        plt.title(r"mass error $\Delta M_t=(M_t-M_0)/M_0$")
