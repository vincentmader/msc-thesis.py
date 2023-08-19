import os
import sys
try:
    sys.path.append(os.path.join("..", "..", "src"))
    from axis import DiscreteMassAxis, DiscreteRadialAxis
    from collision import collision_rate
    from config import Config
    from disk import Disk, DiskRegion
    from dust import particle_radius_from_mass
    from visualization.kernel_v3_2023_08_14.pcolor_matrix_subplot import PcolorMatrixSubplot
    from visualization.kernel_v3_2023_08_14.gridspec_plot import GridspecPlot
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
    xlabel="particle radius $a_j$",
    ylabel="particle radius $a_i$",
    # scales=("log", "log", "lin"),
)

subplots = [s1]

p = GridspecPlot(subplots)
p.render()
