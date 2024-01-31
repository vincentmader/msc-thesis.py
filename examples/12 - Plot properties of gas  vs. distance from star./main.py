import os, sys
from pathlib import Path
import matplotlib.pyplot as plt
try:
    sys.path.append(os.path.join("..", "..", "src"))
    from config import Config, PATH_TO_DARKMODE, PATH_TO_FIGURES
    from constants import AU
    from functions.utils import physics
    from models.disk import Disk, DiskRegion
    from models.axis import DiscreteMassAxis, DiscreteRadialAxis
except ModuleNotFoundError as e:
    raise e

cfg     = Config()
mg      = DiscreteMassAxis(cfg)
rg      = DiscreteRadialAxis(cfg)
rc      = rg.bin_centers
rb      = rg.bin_boundaries
disk    = Disk(cfg, rg, mg)
Sigma_g = disk.gas_surface_density
M_star  = cfg.stellar_mass
distance_to_star = cfg.distance_to_star

T_mid   = disk.midplane_temperature
c_s     = disk.sound_speed
rho_g   = disk.midplane_gas_volume_density
P       = disk.midplane_gas_pressure
u_th    = disk.thermal_velocity
del_ln_P_g_del_ln_r = disk.del_ln_P_g_del_ln_r

disk_region = DiskRegion(cfg, disk)
delr_Sigma_g_nu_g_sqrt_r = disk_region.delr_Sigma_g_nu_g_sqrt_r
u_g = physics.u_g(
    rc, Sigma_g, delr_Sigma_g_nu_g_sqrt_r
)

FIGSIZE = (11, 4)
if cfg.mpl_dark_mode:
    plt.style.use(PATH_TO_DARKMODE)


def plot_1(show_legend=False):
    plt.ylabel(r"gas surface density $\Sigma_g$ [kg/m$^2$]")

    label = r"$\Sigma_g\sim\frac{1}{r}$"
    plt.loglog(rc / AU, Sigma_g, label=label)

    if show_legend:
        plt.legend()


def plot_2(show_legend=False):
    plt.ylabel("midplane temperature $T_{mid}(r)$ [K]")

    label = r"$T_{mid}=(\frac{f}{2}\cdot L_{star}\cdot4\pi\cdot r^2\cdot\sigma_{SB})^{-1/4}$"
    plt.semilogx(rc / AU, T_mid, label=label)

    r = distance_to_star
    L_star = cfg.stellar_luminosity
    phi_fl = cfg.flaring_angle
    T_r = physics.midplane_temperature(r, L_star, phi_fl)
    plt.scatter([r/AU], [T_r], marker="x", color="k", label="position in disk")

    plt.grid(True)
    plt.gca().set_ylim(bottom=0)

    if show_legend:
        plt.legend()


def plot_3(show_legend=False):
    plt.ylabel("sound speed $c_s$ [m/s]")

    label = r"$c_s=\sqrt{\frac{k_BT}{2.3\cdot m_p}}$"
    plt.semilogx(rc / AU, c_s, label=label)

    if show_legend:
        plt.legend()


def plot_4(show_legend=False):
    plt.ylabel(r"gas volume density $\rho_g$ [kg/m$^3$]")

    label = r"gas volume density $\rho_g =\frac{\Sigma_g}{\sqrt{2\pi}\cdot H_p}$"
    plt.loglog(rc / AU, rho_g, label=label)

    if show_legend:
        plt.legend()


def plot_5(show_legend=False):
    plt.ylabel("gas pressure $P$ [Pa]")

    label = r"$P=\rho_g^{mid}\cdot c_s^2$"
    plt.loglog(rc / AU, P, label=label)

    if show_legend:
        plt.legend()


def plot_6(show_legend=False):
    plt.ylabel(r"gas pressure gradient $\frac{\partial\log P_g}{\partial\log r}$")

    label = r"$\frac{\partial\log P_g}{\partial\log r}$"
    plt.semilogx(rc[:-1] / AU, del_ln_P_g_del_ln_r, label=label) # TODO is `[:-1]` the correct slice?

    if show_legend:
        plt.legend()


def plot_7(show_legend=False):
    plt.ylabel(r"radial gas velocity $u_g$ [m/s]")

    label = "TODO"
    plt.semilogx(rc / AU, u_g, label=label) 

    if show_legend:
        plt.legend()


def plot_8(show_legend=False):
    plt.ylabel(r"thermal gas velocity $u_{th}$ [m/s]")

    label = "TODO"
    plt.semilogx(rc / AU, u_th, label=label) 

    if show_legend:
        plt.legend()


PLOTS = {
    "gas_surface_density":  plot_1,
    "midplane_temperature": plot_2,
    "sound_speed":          plot_3,
    "gas_volume_density":   plot_4,
    "gas_pressure":         plot_5,
    "pressure_gradient":    plot_6,
    "radial_gas_velocity":  plot_7,
    "thermal_gas_velocity": plot_8,
}


if __name__ == "__main__":
    os.makedirs(Path(PATH_TO_FIGURES, "12"), exist_ok=True)

    for title, plot in PLOTS.items():
        _, ax = plt.subplots(figsize=FIGSIZE)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

        plt.xlabel("distance from star $r$ [AU]")
        plt.xlim(rc[0] / AU, rc[-1] / AU)
        plt.grid()

        plot()

        path = os.path.join(PATH_TO_FIGURES, "12", f"{title}.pdf")
        plt.savefig(path)
        plt.show()
        plt.close()
