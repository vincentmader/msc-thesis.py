import matplotlib.pyplot as plt

from kernel.mass_conservation import test_mass_conservation
from matplotlib import colors


class KernelMassConservationPlot():

    def __init__(self, cfg, mg, K):
        self.fig = plt.figure()
        self.ax = plt.gca()
        self.K = K
        self.sum_ij = test_mass_conservation(cfg, mg, K)

        self.ax.format_coord = self.custom_format_coord

        COLORBAR_GEOMETRY = [0.5, 0.95, 0.3, 0.01]
        rect = COLORBAR_GEOMETRY

        p = self.draw()
        self.setup_colorbar(p, rect, "horizontal")
        # self.ax.format_coord = format_coord

    def draw(self):
        self.ax.set_aspect('equal', adjustable='box')

        vmin, vmax = 1e-14, 1e14  # TODO Calculate dynamically.
        plt.set_cmap("Reds")
        cmap_norm = colors.LogNorm(vmin=vmin, vmax=vmax)

        return plt.pcolor(
            self.sum_ij,
            norm=cmap_norm
        )

    def show(self):
        plt.show()
        plt.close()

    def custom_format_coord(self, x, y):
        i, j = int(x), int(y)  # TODO Convention?
        sum_ij = self.sum_ij[i, j]
        text = ""
        text += f"sum_k K_kij = {sum_ij:.2}, "
        text += f"{i = }, "
        text += f"{j = } "
        return text

    def setup_colorbar(self, p, rect, orientation):
        # rect = center_rect(rect, orientation)
        ax = self.fig.add_axes(rect)
        self.fig.colorbar(p, cax=ax, orientation=orientation)
