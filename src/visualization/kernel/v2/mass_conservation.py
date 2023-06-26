import matplotlib.pyplot as plt

from kernel.mass_conservation import test_mass_conservation
from matplotlib import colors


class KernelMassConservationPlot():

    def __init__(self, kernel):
        self.fig = plt.figure()
        self.ax = plt.gca()
        self.kernel = kernel

        COLORBAR_GEOMETRY = [0.5, 0.95, 0.3, 0.01]
        rect = COLORBAR_GEOMETRY

        p = self.draw()
        self.setup_colorbar(p, rect, "horizontal")
        # self.ax.format_coord = format_coord

    def draw(self):
        kernel = self.kernel

        self.ax.set_aspect('equal', adjustable='box')

        vmin, vmax = 1e-14, 1e14
        plt.set_cmap("Reds")
        cmap_norm = colors.LogNorm(vmin=vmin, vmax=vmax)

        sum_ij = test_mass_conservation(kernel)
        return plt.pcolor(
            sum_ij,
            norm=cmap_norm
        )

    def show(self):
        plt.show()
        plt.close()

    def setup_colorbar(self, p, rect, orientation):
        # rect = center_rect(rect, orientation)
        ax = self.fig.add_axes(rect)
        self.fig.colorbar(p, cax=ax, orientation=orientation)
