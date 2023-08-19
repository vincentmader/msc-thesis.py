import os
import sys
try:
    sys.path.append(os.path.join("..", "..", "src"))
    from axis import DiscreteMassAxis
    from collision import collision_cross_section
    from config import Config
    from dust import particle_radius_from_mass
    from visualization.vX_2023_08_14_kernel.pcolor_matrix_subplot import PcolorMatrixSubplot
    from visualization.vX_2023_08_14_kernel.gridspec_plot import GridspecPlot
except ModuleNotFoundError as e:
    raise e


cfg = Config()

mg = DiscreteMassAxis(cfg)
mc = mg.grid_cell_centers
rho_s = cfg.dust_particle_density
ac = particle_radius_from_mass(mc, rho_s)

R_coll = collision_cross_section(cfg, mg)

s1 = PcolorMatrixSubplot(
    ac, ac, R_coll,
    title="collision cross section $\sigma_{ij}$",
    xlabel="particle radius $a_j$",
    ylabel="particle radius $a_i$",
    # scales=("log", "log", "lin"),
)

subplots = [s1]

p = GridspecPlot(subplots)
p.render()
