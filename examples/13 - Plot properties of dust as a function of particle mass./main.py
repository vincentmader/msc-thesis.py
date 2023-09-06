import os
import sys
from pathlib import Path

import matplotlib.pyplot as plt
try:
    sys.path.append(os.path.join("..", "..", "src"))
    from axis import DiscreteMassAxis, DiscreteRadialAxis
    from config import Config, PATH_TO_DARKMODE, PATH_TO_FIGURES
    from disk import Disk, DiskRegion
    from dust import particle_radius_from_mass
except ModuleNotFoundError as e:
    raise e

FIGSIZE = (8, 4)

cfg = Config()
mg = DiscreteMassAxis(cfg)
rg = DiscreteRadialAxis(cfg)
disk = Disk(cfg, rg, mg)
disk_region = DiskRegion(cfg, disk)

if cfg.mpl_dark_mode:
    plt.style.use(PATH_TO_DARKMODE)


def plot_1(m, St):
    plt.figure(figsize=FIGSIZE)
    # TODO Include both cases.
    plt.loglog(
        m, St, label=r"$St=\frac{\pi}{2}\cdot\frac{\rho_s\cdot a}{\Sigma_g}$")
    plt.title("Stokes' nr. vs. mass $m$")
    plt.xlabel("mass $m$ [kg]")
    plt.ylabel(r"$St$")
    plt.legend()
    plt.grid()
    path = os.path.join(PATH_TO_FIGURES, "12", "stokes_number.pdf")
    plt.savefig(path)
    plt.show()
    plt.close()


def plot_2(m, Re):
    plt.figure(figsize=FIGSIZE)
    plt.loglog(m, Re, label=r"$Re=\frac{2au}{\nu_{mol}}$")  # TODO Specify 'u'.
    plt.title("Reynold's nr. vs. mass $m$")
    plt.xlabel("mass $m$ [kg]")
    plt.ylabel(r"$Re$")
    plt.legend()
    plt.grid()
    path = os.path.join(PATH_TO_FIGURES, "12", "reynolds_number.pdf")
    plt.savefig(path)
    plt.show()
    plt.close()


def plot_3(m, t_stop):
    plt.figure(figsize=FIGSIZE)
    plt.loglog(
        m, t_stop, label=r"$\tau_{stop}=\frac{\rho_s\cdot a}{\rho_g\cdot u_{th}}$")
    plt.title(r"stopping time $\tau_{stop}$ (Epstein regime)")
    plt.ylabel(r"$\tau_{stop}$ [s]")
    plt.xlabel(r"$m$ [kg]")
    plt.legend()
    plt.grid()
    path = os.path.join(PATH_TO_FIGURES, "12", "stopping_time.pdf")
    plt.savefig(path)
    plt.show()
    plt.close()


if __name__ == "__main__":
    mc = mg.grid_cell_centers
    rho_s = cfg.dust_particle_density
    radii = particle_radius_from_mass(mc, rho_s)
    stopping_times = disk_region.stopping_time(radii, rho_s)
    stokes_nrs = disk_region.stokes_nr(mc, stopping_times, rho_s)
    u = disk_region.v_K  # TODO Is this correct?
    reynolds_nrs = disk_region.reynolds_nr(mc, u)

    os.makedirs(Path(PATH_TO_FIGURES, "12"), exist_ok=True)
    plot_1(mc, stokes_nrs)
    plot_2(mc, reynolds_nrs)
    plot_3(mc, stopping_times)