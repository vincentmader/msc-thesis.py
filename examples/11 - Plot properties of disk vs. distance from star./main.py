import os, sys
from pathlib import Path
import matplotlib.pyplot as plt
try:
    sys.path.append(os.path.join("..", "..", "src"))
    from config import Config, PATH_TO_DARKMODE, PATH_TO_FIGURES
    from constants import AU
    from models.disk import Disk
    from models.axis import DiscreteMassAxis, DiscreteRadialAxis
    from functions.utils.physics import kepler_frequency
except ModuleNotFoundError as e:
    raise e

cfg = Config()
mg = DiscreteMassAxis(cfg)
rg = DiscreteRadialAxis(cfg)
rc = rg.bin_centers
rb = rg.bin_boundaries
disk = Disk(cfg, rg, mg)
Sigma_g = disk.gas_surface_density
M_star = cfg.stellar_mass

FIGSIZE = (10, 4)
if cfg.mpl_dark_mode:
    plt.style.use(PATH_TO_DARKMODE)


def plot_1(show_legend=False):
    H_p = disk.scale_height
    plt.title("disk scale height $H_p$")
    plt.ylabel("$H_p$ [AU]")
    plt.loglog(rc / AU, H_p / AU, label=r"$H_p=\frac{c_s}{\Omega_K}$")
    plt.legend()


def plot_2(show_legend=False):
    Omega_K = kepler_frequency(rc, M_star)
    plt.title("Kepler frequency $\Omega_K$")
    plt.ylabel("$\Omega_K$ [1/s]")
    label = r"$\Omega_K=\sqrt{\frac{G\cdot M_{star}}{r^3}}$"
    plt.loglog(rc / AU, Omega_K, label=label)
    plt.legend()


PLOTS = {
    "disk_scale_height": plot_1,
    "kepler_frequency": plot_2,
}


if __name__ == "__main__":
    os.makedirs(Path(PATH_TO_FIGURES, "02"), exist_ok=True)

    for title, plot in PLOTS.items():
        _, ax = plt.subplots(figsize=FIGSIZE)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

        plt.xlabel("distance from star $r$ [AU]")
        plt.xlim(rc[0] / AU, rc[-1] / AU)
        plt.grid()

        plot()

        path = os.path.join(PATH_TO_FIGURES, "02", f"{title}.pdf")
        plt.savefig(path)
        plt.show()
        plt.close()
