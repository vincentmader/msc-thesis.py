import os, sys
from pathlib import Path
import matplotlib.pyplot as plt
try:
    sys.path.append(os.path.join("..", "..", "src"))
    from config import Config, PATH_TO_DARKMODE, PATH_TO_FIGURES
    from constants import AU
    from functions.utils import physics
    from models.disk import Disk
    from models.axis import DiscreteMassAxis, DiscreteRadialAxis
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
distance_to_star = cfg.distance_to_star

FIGSIZE = (10, 4)
if cfg.mpl_dark_mode:
    plt.style.use(PATH_TO_DARKMODE)


def assure_existence_of_figure_directory():
    path_to_figures = Path(PATH_TO_FIGURES, "12")
    os.makedirs(path_to_figures, exist_ok=True)


def plot_1():
    plt.title("gas surface density $\Sigma_g$")
    plt.ylabel(r"$\Sigma_g$ [kg/m$^2$]")
    plt.xlabel("distance from star $r$ [AU]")
    label = r"$\Sigma_g\sim\frac{1}{r}$"
    plt.loglog(rc / AU, Sigma_g, label=label)
    plt.legend()


def plot_2():
    T_mid = disk.midplane_temperature
    plt.ylabel("midplane temperature $T_{mid}(r)$ [K]")
    plt.xlabel("distance from star $r$ [AU]")

    label= "midplane temperature"
    # label = r"$T_{mid}=(\frac{f}{2}\cdot L_{star}\cdot4\pi\cdot r^2\cdot\sigma_{SB})^{-1/4}$"
    plt.semilogx(rc / AU, T_mid, label=label)

    r = distance_to_star
    L_star = cfg.stellar_luminosity
    phi_fl = cfg.flaring_angle
    T_r = physics.midplane_temperature(r, L_star, phi_fl)
    plt.scatter([r/AU], [T_r], marker="x", color="k", label="position in disk")

    plt.legend()
    plt.grid(True)
    plt.gca().set_ylim(bottom=0)


def plot_3():
    c_s = disk.sound_speed
    plt.title("sound speed $c_s$")
    plt.ylabel("$c_s$ [m/s]")
    plt.xlabel("distance from star $r$ [AU]")
    label = r"$c_s=\sqrt{\frac{k_BT}{2.3\cdot m_p}}$"
    plt.loglog(rc / AU, c_s, label=label)
    plt.legend()


def plot_4():
    rho_g = disk.midplane_gas_volume_density
    plt.title(r"gas volume density $\rho_g$")
    plt.ylabel(r"$\rho_g$ [kg/m$^3$]")
    plt.xlabel("distance from star $r$ [AU]")
    label = r"$\rho_g =\frac{\Sigma_g}{\sqrt{2\pi}\cdot H_p}$"
    plt.loglog(rc / AU, rho_g, label=label)
    plt.legend()


def plot_5():
    P = disk.midplane_gas_pressure
    plt.title("midplane gas pressure $P$")
    plt.ylabel("$P$ [Pa]")
    plt.xlabel("distance from star $r$ [AU]")
    plt.loglog(rc / AU, P, label=r"$P=\rho_g^{mid}\cdot c_s^2$")
    plt.legend()


def plot_6():
    del_ln_P_g_del_ln_r = disk.del_ln_P_g_del_ln_r
    plt.title("Logarithmic pressure gradient")
    plt.ylabel(r"$\frac{\partial\log P_g}{\partial\log r}$")
    plt.xlabel("distance from star $r$ [AU]")
    label = r"$\frac{\partial\log P_g}{\partial\log r}$"
    plt.semilogx(rc[:-1] / AU, del_ln_P_g_del_ln_r, label=label) # TODO is `[:-1]` the correct slice?
    plt.legend()


def plot_7():
    from models.disk import DiskRegion
    disk_region = DiskRegion(cfg, disk)
    delr_Sigma_g_nu_g_sqrt_r = disk_region.delr_Sigma_g_nu_g_sqrt_r
    v_r = physics.u_g(
        rc, Sigma_g, delr_Sigma_g_nu_g_sqrt_r
    )

    plt.title("Radial gas velocity [m/s]")
    plt.xlabel("distance from star $r$ [AU]")
    plt.ylabel(r"$u_g$")
    plt.semilogx(rc / AU, v_r) 


def plot_8():
    T_mid = disk.midplane_temperature
    c_s = physics.sound_speed(rc, T_mid)

    plt.title("Thermal gas velocity [m/s]")
    plt.xlabel("distance from star $r$ [AU]")
    plt.ylabel(r"$c_s$")
    plt.semilogx(rc / AU, c_s) 


def create_figure(plotter_function):
    _, ax = plt.subplots(figsize=FIGSIZE)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    plt.grid()
    plotter_function()
    plt.xlim(rc[0] / AU, rc[-1] / AU)

    path = os.path.join(PATH_TO_FIGURES, "12", f"{title}.pdf")
    plt.savefig(path)
    plt.show()
    plt.close()


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
    assure_existence_of_figure_directory()
    for title, plotter_function in PLOTS.items():
        create_figure(plotter_function)
