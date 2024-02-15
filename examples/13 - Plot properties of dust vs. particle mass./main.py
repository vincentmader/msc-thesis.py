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

FIGSIZE = (11, 4)

cfg = Config()
mg = DiscreteMassAxis(cfg)
rg = DiscreteRadialAxis(cfg)
disk = Disk(cfg, rg, mg)
disk_region = DiskRegion(cfg, disk)

rho_g_mid   = disk_region.midplane_gas_volume_density
u_th        = disk_region.thermal_velocity
mc          = disk_region.mg.bin_centers
u           = disk_region.kepler_velocity  # TODO Is this correct?
Sigma_g     = disk_region.gas_surface_density
nu_mol      = disk_region.viscosity
rho_s       = cfg.dust_particle_density
a           = particle_radius_from_mass(mc, rho_s)

Re          = physics.reynolds_nr(a, u, nu_mol)
t_stop      = physics.stopping_time(rho_s, a, rho_g_mid, u_th)
St          = physics.stokes_nr(rho_s, a, Sigma_g)

if cfg.mpl_dark_mode:
    plt.style.use(PATH_TO_DARKMODE)


def plot_1(show_legend=False):
    plt.ylabel(r"Stokes Number $St$")

    label = r"$St=\frac{\pi}{2}\cdot\frac{\rho_s\cdot a}{\Sigma_g}$"
    plt.loglog(mc, St, label=label)

    if show_legend:
        plt.legend()


def plot_2(show_legend=False):
    plt.ylabel(r"Reynolds Number $Re$")

    label = r"$Re=\frac{2au}{\nu_{mol}}$"
    plt.loglog(mc, Re, label=label)  # TODO Specify 'u'.

    if show_legend:
        plt.legend()


def plot_3(show_legend=False):
    plt.ylabel(r"Stopping Time $\tau_{stop}$ [s]")

    label = r"$\tau_{stop}=\frac{\rho_s\cdot a}{\rho_g\cdot u_{th}}$ (Epstein regime)"
    plt.loglog(mc, t_stop, label=label)

    if show_legend:
        plt.legend()


if __name__ == "__main__":
    os.makedirs(Path(PATH_TO_FIGURES, "13"), exist_ok=True)

    PLOTS = {
        "stokes_nr":        plot_1,
        "reynolds_nr":      plot_2,
        "stopping_time":    plot_3,
    }

    for plot_id, plot in PLOTS.items():
        _, ax = plt.subplots(figsize=FIGSIZE)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

        plt.xlabel(r"dust particle mass $m$ [kg]")
        plt.xlim(mc[0], mc[-1])
        plt.grid()

        plot()

        path = os.path.join(PATH_TO_FIGURES, "13", f"{plot_id}.pdf")
        plt.savefig(path)
        plt.show()
        plt.close()
