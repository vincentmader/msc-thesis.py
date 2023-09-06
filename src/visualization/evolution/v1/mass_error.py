import matplotlib.pyplot as plt

from disk.disk import disk_mass_error


FIGSIZE = (6, 4)


class DiskMassErrorPlot():

    def __init__(
        self,
        time,
        disk_masses,
        figsize=None,
    ):
        self.x = time
        self.y = disk_mass_error(disk_masses, disk_masses[0])

        if figsize is None:
            figsize = FIGSIZE
        self.fig = plt.figure(figsize=figsize)
        self.ax_1 = plt.gca()

    def draw(self):
        x = self.x / (365 * 24 * 60 * 60)
        y = self.y

        y = y[:-1] / (x[1:] - x[:-1])
        x = x[:-1]

        self.ax_1.loglog(x, y, 'r', label="pos.")
        self.ax_1.loglog(x, -y, 'b', label="neg.")
        # TODO ^ Draw as semilogy for linear time axis.
        plt.title(r"mass error $\frac{m-m_0}{m_0}$")
        plt.ylabel("$n_i$ [kg$^{-1}$ $m^{-3}$ $s^{-1}}$]")
        plt.xlabel("time [years]")
        plt.ylabel(r"$\frac{m-m_0}{m_0}$")
        plt.legend()
