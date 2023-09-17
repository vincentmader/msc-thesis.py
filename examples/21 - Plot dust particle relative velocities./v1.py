import os, sys
import matplotlib.pyplot as plt
from tqdm import tqdm
try:
    sys.path.append(os.path.join("..", "..", "src"))
    from axis import DiscreteMassAxis, DiscreteRadialAxis
    from config import Config, PATH_TO_DARKMODE
    from disk import Disk, DiskRegion
    from dust import particle_radius_from_mass
    from dust.relative_velocity import dv_azimuthal
    from dust.relative_velocity import dv_brownian_motion
    from dust.relative_velocity import dv_differential_settling
    from dust.relative_velocity import dv_radial_drift
    from dust.relative_velocity import dv_turbulence
    from dust.relative_velocity import relative_velocity
except ModuleNotFoundError as e:
    raise e

# Define configuration.
cfg = Config(
    mass_max_value=1e6,
    mass_resolution=50,
    # relative_velocity_components=[
    #     "brownian_motion",
    #     "differential_settling",
    #     "radial_drift",
    #     "turbulence",
    #     "azimuthal",
    # ],
)

# Setup pyplot figure.
FIGSIZE_MULTI = (13, 7)
FIGSIZE_SINGLE = (8, 6)
SLIDER_POSITION = [0.05, 0.3, 0.02, 0.4]
if cfg.mpl_dark_mode:
    plt.style.use(PATH_TO_DARKMODE)

# Define discrete axis for radial distance from star, as well as for mass.
rg = DiscreteRadialAxis(cfg)
mg = DiscreteMassAxis(cfg)
N_m = mg.N
masses = mg.bin_centers  # TODO Use bounds or centers?

rho_s = cfg.dust_particle_density

# Calculate particle radii from masses.
radii = particle_radius_from_mass(masses, rho_s)

# Define disk, the position of interest in it, & the disk properties there.
disk = Disk(cfg, rg, mg)
disk_region = DiskRegion(cfg, disk)

# Calculate relative velocities.
dv_br = dv_brownian_motion(cfg, disk, disk_region)
dv_ds = dv_differential_settling(cfg, disk, disk_region)
dv_rd = dv_radial_drift(cfg, disk, disk_region)
dv_az = dv_azimuthal(cfg, disk, disk_region)
dv_tu = dv_turbulence(cfg, disk, disk_region)
dv = relative_velocity(cfg, disk, disk_region)

# Define plot setups.
plot_setups = [
    (dv_br, "BR"),
    (dv_az, "AZ"),
    (dv_rd, "RD"),
    (dv_tu, "TU"),
    (dv_ds, "DS"),
    (dv, "tot"),
]


# def upper_half(m):
#     res = m
#     for i in range(len(m)):
#         for j in range(len(m)):
#             if i < j:
#                 res[i, j] = None
#     return res


# def foo():
#     ticks, locs = [], []
#     ns = []  # n = order of magnitude
#     for i, m in enumerate(mg.bin_boundaries):
#         n = int(np.log10(m))
#         if n not in ns:
#             ns.append((i, n))
#     for i, n in ns:
#         m = 10**n
#         locs.append(i)
#         ticks.append(m)
#     ticks = ['%.0e' % i for i in ticks]
#     return locs[::10], ticks[::10]


def plot(dv, title):
    ax = plt.gca()
    plt.title(r"$\Delta v_{" + title + r"}$")
    plt.set_cmap("Reds")
    # plt.ylabel("$i$", rotation=0)
    # plt.xlabel("$j$")

    # dv = upper_half(dv)

    # ticks = radii
    # ticks = ticks[::10]
    # ticks = [np.log(i) for i in ticks]
    # ticks = [np.round(i, 1) for i in ticks]
    # locs = range(N_m)[::10]
    # plt.xticks(locs, ticks)
    # plt.yticks(locs, ticks)

    # plt.pcolor(np.log(dv))
    plt.pcolor(dv)
    plt.colorbar()
    plt.axis("scaled")
    # TODO Add axis labels.

    def format_num(x):
        return "%.2e" % x

    def format_coord(x, y):
        i = int(x)
        j = int(y)
        m_i = masses[i]
        m_j = masses[j]
        a_i = particle_radius_from_mass(m_i, rho_s)
        a_j = particle_radius_from_mass(m_j, rho_s)
        t_s_i = disk_region.stopping_time(a_i, rho_s)
        t_s_j = disk_region.stopping_time(a_j, rho_s)
        St_i = disk_region.stokes_nr(m_i, t_s_i, rho_s)
        St_j = disk_region.stokes_nr(m_j, t_s_j, rho_s)
        m_i, m_j = format_num(m_i), format_num(m_j)
        a_i, a_j = format_num(a_i), format_num(a_j)
        St_i, St_j = format_num(St_i), format_num(St_j)
        dv_ij = dv[i, j]
        dv_ij = format_num(dv_ij)
        text = f"{i=},  \tm_i={m_i} kg,  \ta_i={a_i} m,  \tSt_i={St_i}\n"
        text += f"{j=},  \tm_j={m_j} kg,  \ta_j={a_j} m,  \tSt_j={St_j},  \tdv_ij={dv_ij} m/s"
        return text

    ax.format_coord = format_coord


def plot_separately():
    for dv, title in tqdm(plot_setups):
        plt.figure(
            figsize=FIGSIZE_SINGLE,
            constrained_layout=True,
        )
        plot(dv, title)
        # path = os.path.join(PATH_TO_FIGURES, "13", f"dv_{title}.pdf")
        # plt.savefig(path)
        plt.close()


def plot_together():
    fig, ax = plt.subplots(figsize=FIGSIZE_MULTI)
    for idx, (dv, title) in tqdm(enumerate(plot_setups)):
        plt.subplot(2, 3, idx + 1)
        plot(dv, title)
    # path = os.path.join(PATH_TO_FIGURES, "13", "relative_velocities.pdf")
    # plt.savefig(path)
    plt.show()
    plt.close()


def main():
    # os.makedirs("../../figures/13", exist_ok=True)
    plot_separately()
    plot_together()
