import os, sys
from pathlib import Path
import matplotlib.pyplot as plt
try:
    sys.path.append(os.path.join("..", "..", "src"))
    from models.axis import DiscreteMassAxis, DiscreteRadialAxis
    from config import Config, PATH_TO_DARKMODE, PATH_TO_FIGURES
    from models.disk import Disk, DiskRegion
    from functions.dust import particle_radius_from_mass
    from functions.utils import physics
except ModuleNotFoundError as e:
    raise e

FIGSIZE = (10, 4)

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
    rho_g_mid   = disk_region.midplane_gas_volume_density
    u_th        = disk_region.thermal_velocity
    mc          = disk_region.mg.bin_centers
    u           = disk_region.kepler_velocity  # TODO Is this correct?
    Sigma_g     = disk_region.gas_surface_density
    nu_mol      = disk_region.viscosity
    rho_s       = cfg.dust_particle_density
    a           = particle_radius_from_mass(mc, rho_s)

    os.makedirs(Path(PATH_TO_FIGURES, "12"), exist_ok=True)

    reynolds_nrs = physics.reynolds_nr(a, u, nu_mol)
    plot_2(mc, reynolds_nrs)

    stopping_times = physics.stopping_time(rho_s, a, rho_g_mid, u_th)
    plot_3(mc, stopping_times)

    stokes_nrs = physics.stokes_nr(rho_s, a, Sigma_g)
    plot_1(mc, stokes_nrs)
