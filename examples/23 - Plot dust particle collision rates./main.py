import os, sys
from pathlib import Path
try:
    sys.path.append(os.path.join("..", "..", "src"))
    from axis import DiscreteMassAxis, DiscreteRadialAxis
    from collision import collision_rate
    from config import Config, PATH_TO_FIGURES
    from disk import Disk, DiskRegion
    from dust import particle_radius_from_mass
    from visualization.base import GridspecPlot, PcolorMatrixSubplot
except ModuleNotFoundError as e:
    raise e


cfg = Config()

rg = DiscreteRadialAxis(cfg)
mg = DiscreteMassAxis(cfg)

mc = mg.grid_cell_centers
rho_s = cfg.dust_particle_density
ac = particle_radius_from_mass(mc, rho_s)

disk = Disk(cfg, rg, mg)
disk_region = DiskRegion(cfg, disk)

R_coll = collision_rate(cfg, disk, disk_region)

s1 = PcolorMatrixSubplot(
    ac, ac, R_coll,
    title="collision rate $R_{coll}$",
    xlabel="particle radius $a_j$ [m]",
    ylabel="particle radius $a_i$ [m]",
    # scales=("log", "log", "lin"),
    cmap="Blues",
)

path_to_figures = Path(PATH_TO_FIGURES, "23")
os.makedirs(path_to_figures, exist_ok=True)
path_to_outfile = Path(path_to_figures, "collision_rate.pdf")

p = GridspecPlot([s1])
p.render(
    save_plot=True,
    path_to_outfile=path_to_outfile
)
