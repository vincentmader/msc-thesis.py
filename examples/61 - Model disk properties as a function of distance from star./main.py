import os
import sys
import matplotlib.pyplot as plt
try:
    sys.path.append(os.path.join("..", "..", "src"))
    from disk import Disk, MassGrid
    from utils.axis import DiscreteAxis
    from utils.plotting import plt_show_then_close
    from utils.physics import kepler_frequency
    from config import Config, PATH_TO_DARKMODE, PATH_TO_FIGURES
    from constants import AU
except ModuleNotFoundError as e:
    raise e

FIGSIZE = (10, 5)
cfg = Config()
mg = MassGrid(cfg)
r_min, r_max = cfg.radial_min_value, cfg.radial_max_value
N_r, scale = cfg.radial_resolution, cfg.radial_axis_scale
rg = DiscreteAxis(r_min, r_max, N_r, scale)
r = rg.grid_cell_centers()
r_bounds = rg.grid_cell_boundaries()
disk = Disk(cfg, rg, mg)
Sigma_g = disk.gas_surface_density(r_bounds)
M_star = cfg.stellar_mass

if cfg.mpl_dark_mode:
    plt.style.use(PATH_TO_DARKMODE)


def plot_1():
    plt.title("gas surface density $\Sigma_g$")
    plt.ylabel("$\Sigma_g$ [kg/m$^2$]")
    plt.xlabel("distance from star $r$ [AU]")
    label = r"$\Sigma_g\sim\frac{1}{r}$"
    plt.loglog(r / AU, Sigma_g, label=label)
    plt.legend()


def plot_2():
    T_mid = disk.midplane_temperature(r)
    plt.title("midplane temperature $T_{mid}$")
    plt.ylabel("$T_{mid}$ [K]")
    plt.xlabel("distance from star $r$ [AU]")
    label = r"$T_{mid}=\frac{f}{2}\cdot L_{star}\cdot(4\pi\cdot r^2\cdot\sigma_{SB})^{-1/4}$"
    plt.loglog(r / AU, T_mid, label=label)
    plt.legend()


def plot_3():
    c_s = disk.sound_speed(r)
    plt.title("sound speed $c_s$")
    plt.ylabel("$c_s$ [m/s]")
    plt.xlabel("distance from star $r$ [AU]")
    label = r"$c_s=\sqrt{\frac{k_BT}{2.3\cdot m_p}}$"
    plt.loglog(r / AU, c_s, label=label)
    plt.legend()


def plot_4():
    H_p = disk.scale_height(r)
    plt.title("disk scale height $H_p$")
    plt.ylabel("$H_p$ [AU]")
    plt.xlabel("distance from star $r$ [AU]")
    plt.loglog(r / AU, H_p / AU, label=r"$H_p=\frac{c_s}{\Omega_K}$")
    plt.legend()


def plot_5():
    rho_g = disk.midplane_gas_volume_density(r, Sigma_g)
    plt.title(r"gas volume density $\rho_g$")
    plt.ylabel(r"$\rho_g$ [kg/m$^3$]")
    plt.xlabel("distance from star $r$ [AU]")
    label = r"$\rho_g =\frac{\Sigma_g}{\sqrt{2\pi}\cdot H_p}$"
    plt.loglog(r / AU, rho_g, label=label)
    plt.legend()


def plot_6():
    P = disk.midplane_gas_pressure(r, Sigma_g)
    plt.title("midplane gas pressure $P$")
    plt.ylabel("$P$ [Pa]")
    plt.xlabel("distance from star $r$ [AU]")
    plt.loglog(r / AU, P, label=r"$P=\rho_g^{mid}\cdot c_s^2$")
    plt.legend()


def plot_7():
    Omega_K = kepler_frequency(r, M_star)
    plt.title("Kepler frequency $\Omega_K$")
    plt.ylabel("$\Omega_K$ [1/s]")
    plt.xlabel("distance from star $r$ [AU]")
    label = r"$\Omega_K=\sqrt{\frac{G\cdot M_{star}}{r^3}}$"
    plt.loglog(r / AU, Omega_K, label=label)
    plt.legend()


def plot_8():
    del_ln_P_g_del_ln_r = disk.del_ln_P_g_del_ln_r(r, Sigma_g)
    plt.title("Logarithmic pressure gradient")
    plt.ylabel(r"$\frac{\partial\log P_g}{\partial\log r}$")
    plt.xlabel("distance from star $r$ [AU]")
    label = r"$\frac{\partial\log P_g}{\partial\log r}$"
    plt.semilogx(r / AU, del_ln_P_g_del_ln_r, label=label)
    plt.legend()

PLOTS = {
    "gas_surface_density": plot_1,
    "midplane_temperature": plot_2,
    "sound_speed": plot_3,
    "disk_scale_height": plot_4,
    "gas_volume_density": plot_5,
    "gas_pressure": plot_6,
    "kepler_frequency": plot_7,
    "pressure_gradient": plot_8,
}

if __name__ == "__main__":
    for title, plot in PLOTS.items():
        fig, ax = plt.subplots(figsize=FIGSIZE)
        plot()
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        path = os.path.join(PATH_TO_FIGURES, "02", f"{title}.pdf")
        plt.savefig(path)
        plt_show_then_close()
