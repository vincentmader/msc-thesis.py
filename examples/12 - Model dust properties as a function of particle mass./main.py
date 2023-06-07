import os
import sys
import matplotlib.pyplot as plt
try:
    sys.path.append(os.path.join("..", "..", "src"))
    from config import Config, PATH_TO_DARKMODE, PATH_TO_FIGURES
    from disk import MassGrid, Disk, DiskRegion, RadialGrid
    from utils.plotting import plt_show_then_close
    from disk.dust_particle import particle_radius_from_mass
except ModuleNotFoundError as e:
    raise e

FIGSIZE = (10, 5)

cfg = Config()
mg = MassGrid(cfg)
rg = RadialGrid(cfg)
disk = Disk(cfg, rg, mg)
disk_region = DiskRegion(cfg, disk)

if cfg.mpl_dark_mode:
    plt.style.use(PATH_TO_DARKMODE)


def plot_1(m, St):
    plt.figure(figsize=FIGSIZE)
    plt.loglog(m, St, label=r"$St=\frac{\pi}{2}\cdot\frac{\rho_s\cdot a}{\Sigma_g}$") # TODO Include both cases.
    plt.title("Stokes' nr. vs. mass $m$")
    plt.xlabel("mass $m$ [kg]")
    plt.ylabel(r"$St$")
    plt.legend()
    path = os.path.join(PATH_TO_FIGURES, "03", "stokes_number.pdf")
    plt.savefig(path)
    plt_show_then_close()

def plot_2(m, Re):
    plt.figure(figsize=FIGSIZE)
    plt.loglog(m, Re, label=r"$Re=\frac{2au}{\nu_{mol}}$") # TODO Specify 'u'.
    plt.title("Reynold's nr. vs. mass $m$")
    plt.xlabel("mass $m$ [kg]")
    plt.ylabel(r"$Re$")
    plt.legend()
    path = os.path.join(PATH_TO_FIGURES, "03", "reynolds_number.pdf")
    plt.savefig(path)
    plt_show_then_close()

def plot_3(m, t_stop):
    plt.figure(figsize=FIGSIZE)
    plt.loglog(m, t_stop, label=r"$\tau_{stop}=\frac{\rho_s\cdot a}{\rho_g\cdot u_{th}}$")
    plt.title(r"stopping time $\tau_{stop}$ (Epstein regime)")
    plt.ylabel(r"$\tau_{stop}$ [s]")
    plt.xlabel(r"$m$ [kg]")
    plt.legend()
    path = os.path.join(PATH_TO_FIGURES, "03", "stopping_time.pdf")
    plt.savefig(path)
    plt_show_then_close()


if __name__ == "__main__":
    masses = mg.grid_cell_boundaries()  # TODO Use boundaries or centers?
    radii = particle_radius_from_mass(masses)
    stopping_times = disk_region.stopping_time(radii)
    stokes_nrs = disk_region.stokes_nr(masses, stopping_times)
    u = disk_region.v_K  # TODO Is this correct?
    reynolds_nrs = disk_region.reynolds_nr(masses, u)

    plot_1(masses, stokes_nrs)
    plot_2(masses, reynolds_nrs)
    plot_3(masses, stopping_times)
