import os, sys
from pathlib import Path
try:
    sys.path.append(os.path.join("..", "..", "src"))
    from axis import DiscreteMassAxis, DiscreteRadialAxis
    from config import Config
    from config import PATH_TO_FIGURES
    from disk import Disk, DiskRegion
    from dust.relative_velocity import dv_azimuthal
    from dust.relative_velocity import dv_brownian_motion
    from dust.relative_velocity import dv_differential_settling
    from dust.relative_velocity import dv_radial_drift
    from dust.relative_velocity import dv_turbulence
    from dust.relative_velocity import relative_velocity
    from visualization.base import GridspecPlot, PcolorMatrixSubplot
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
dv_ds = dv_differential_settling(cfg, disk, disk_region)
dv_rd = dv_radial_drift(cfg, disk, disk_region)
dv_az = dv_azimuthal(cfg, disk, disk_region)
dv_tu = dv_turbulence(cfg, disk, disk_region)
dv = relative_velocity(cfg, disk, disk_region)

plot_setups = [
    (dv_br, "BR"),
    (dv_az, "AZ"),
    (dv_rd, "RD"),
    (dv_tu, "TU"),
    # (dv_ds, "DS"),
    (dv, "tot"),
]


if __name__ == "__main__":
    for dv, dv_id in plot_setups:
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
        )
    
        p = GridspecPlot([s1])
        p.render(
            save_plot=True,
            path_to_outfile=path_to_outfile,
        )
