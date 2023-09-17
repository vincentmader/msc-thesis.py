import os
import sys
from pathlib import Path
try:
    sys.path.append(os.path.join("..", "..", "src"))
    from axis import DiscreteMassAxis
    from collision import collision_cross_section
    from config import Config, PATH_TO_FIGURES
    from dust import particle_radius_from_mass
    from visualization.base import GridspecPlot, PcolorMatrixSubplot
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
    xlabel="particle radius $a_j$ [m]",
    ylabel="particle radius $a_i$ [m]",
    # scales=("log", "log", "lin"),
    cmap="Blues",
)

path_to_figures = Path(PATH_TO_FIGURES, "22")
os.makedirs(path_to_figures, exist_ok=True)
path_to_outfile = Path(path_to_figures, "collision_cross_section.pdf")

p = GridspecPlot([s1])
p.render(
    save_plot=True,
    path_to_outfile=path_to_outfile
)
