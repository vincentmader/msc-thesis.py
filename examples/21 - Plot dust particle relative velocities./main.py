import os, sys
from pathlib import Path
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, Normalize
import numpy as np
try:
    sys.path.append(os.path.join("..", "..", "src"))
    from models.axis import DiscreteMassAxis, DiscreteRadialAxis
    from config import Config
    from config import PATH_TO_FIGURES
    from models.disk import Disk, DiskRegion
    from functions.dust.relative_velocity import dv_azimuthal
    from functions.dust.relative_velocity import dv_brownian_motion
    from functions.dust.relative_velocity import dv_differential_settling
    from functions.dust.relative_velocity import dv_radial_drift
    from functions.dust.relative_velocity import dv_turbulence
    from functions.dust.relative_velocity import relative_velocity
    from models.plotting.base import GridspecPlot, PcolorMatrixSubplot
except ModuleNotFoundError as e:
    raise e

cfg = Config(
    mass_resolution=200,
    mass_max_value=1e12,
)

rg = DiscreteRadialAxis(cfg)
mg = DiscreteMassAxis(cfg)
mc = mg.bin_centers
ac = mg.particle_radii

# Define disk, the position of interest in it, & the disk properties there.
disk = Disk(cfg, rg, mg)
disk_region = DiskRegion(cfg, disk)

# Calculate relative velocities.
dv_br = dv_brownian_motion(cfg, disk, disk_region)
# dv_ds = dv_differential_settling(cfg, disk, disk_region)
dv_rd = dv_radial_drift(cfg, disk, disk_region)
dv_az = dv_azimuthal(cfg, disk, disk_region)
dv_tu = dv_turbulence(cfg, disk, disk_region)
dv = relative_velocity(cfg, disk, disk_region)

plot_setups = [
    (dv_br, "BR",  (0, 1.4e-3)),
    (dv_az, "AZ",  (0, 25)),
    (dv_rd, "RD",  (0, 25)),
    (dv_tu, "TU",  (0, 25)),
    # (dv_ds, "DS"),
    (dv,    "tot", (0, 40)),
]


def v1():
    for dv, dv_id, z_limits in plot_setups:
        path_to_outfiles = Path(PATH_TO_FIGURES, "21")
        os.makedirs(path_to_outfiles, exist_ok=True)
        filename = f"dv_{dv_id}.pdf"
        path_to_outfile = Path(path_to_outfiles, filename)

        title = r"relative velocity $\Delta v_{ij}^{" + dv_id + "}$ [m/s]"
        s1 = PcolorMatrixSubplot(
            ac, ac, dv,
            title=title,
            xlabel="particle radius $a_j$ [m]",
            ylabel="particle radius $a_i$ [m]",
            axis_scales=("log", "log", "lin"),
            cmap="Reds",
            z_limits=z_limits,
        )
    
        p = GridspecPlot([s1], 
            # figsize=(10, 12) if dv_id == "tot" else None,
        )
        p.render(
            save_plot=True,
            path_to_outfile=path_to_outfile,
        )

def v2():
    path_to_outfiles = Path(PATH_TO_FIGURES, "21")
    os.makedirs(path_to_outfiles, exist_ok=True)
    filename = f"dv_tot_lin+log.pdf"
    path_to_outfile = Path(path_to_outfiles, filename)

    s1 = PcolorMatrixSubplot(
        ac, ac, dv,
        title=r"$\Delta v_{ij}^{tot}$ [m/s]",
        xlabel="particle radius $a_j$ [m]",
        ylabel="particle radius $a_i$ [m]",
        axis_scales=("log", "log", "lin"),
        cmap="Reds",
        z_limits=(0, 40),
        tight_layout=True,
    )
    s2 = PcolorMatrixSubplot(
        ac, ac, dv,
        title=r"$\log\left(\Delta v_{ij}^{tot}\ /\ \frac{m}{s}\right)$",
        xlabel="particle radius $a_j$ [m]",
        ylabel="particle radius $a_i$ [m]",
        axis_scales=("log", "log", "log"),
        cmap="Reds",
        z_limits=(5e-3, 4e1),
        tight_layout=True,
    )

    p = GridspecPlot([s1, s2],
        # figsize=(10, 12) if dv_id == "tot" else None,
    )
    p.render(
        save_plot=True,
        path_to_outfile=path_to_outfile,
    )


def v3():
    path_to_outfiles = Path(PATH_TO_FIGURES, "21")
    os.makedirs(path_to_outfiles, exist_ok=True)
    filename = f"dv_all_2x2.pdf"
    path_to_outfile = Path(path_to_outfiles, filename)

    subplots = []
    for idx, (dv, dv_id, z_limits) in enumerate(plot_setups[:-1]):

        title = r"$\Delta v_{ij}^{" + dv_id + "}$ [m/s]"
        s = PcolorMatrixSubplot(
            ac, ac, dv,
            title=title,
            xlabel="particle radius $a_j$ [m]" if idx > len(plot_setups)-4 else "",
            ylabel="particle radius $a_i$ [m]" if idx % 2 == 0 else "",
            axis_scales=("log", "log", "lin"),
            cmap="Reds",
            z_limits=z_limits,
            tight_layout=True,
        )
        subplots.append(s)
    
    p = GridspecPlot(subplots,
        figsize=(10, 10),
        gridspec_dimensions=(2, 2),
    )
    p.render(
        save_plot=True,
        path_to_outfile=path_to_outfile,
    )

def v4():

    plt.figure(figsize=(10, 8))

    for idx, (dv, dv_id, z_limits) in enumerate(plot_setups[:-1]):
        ax = plt.subplot(2, 2, idx + 1)

        plt.title(r"$\Delta v_{ij}^{" + dv_id + "}$ [m/s]")
        im = plt.pcolormesh(
            ac, ac, dv, 
            cmap="Reds", rasterized=True,
            norm=Normalize(vmin=z_limits[0], vmax=z_limits[1])
        )
        plt.axis("scaled")
        ax.set_xscale("log")
        ax.set_yscale("log")
        # if idx % 2 == 0:
        plt.ylabel("particle radius $a_i$")
        # if idx + 3 > len(plot_setups[:-1]):
        plt.xlabel("particle radius $a_j$")
        plt.colorbar(im, fraction=0.046, pad=0.04)

    plt.tight_layout()
    path_to_outfiles = Path(PATH_TO_FIGURES, "21")
    os.makedirs(path_to_outfiles, exist_ok=True)
    filename = f"dv_all_2x2.pdf"
    plt.savefig(Path(path_to_outfiles, filename))

    plt.show()
    plt.close()

if __name__ == "__main__":
    # v1()
    # v2()
    # v3()
    v4()
