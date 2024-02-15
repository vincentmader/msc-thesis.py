import os, sys
from pathlib import Path
try:
    sys.path.append(os.path.join("..", "..", "src"))
    from models.axis import DiscreteMassAxis
    from functions.dust.collision import collision_cross_section
    from config import Config, PATH_TO_FIGURES
    from models.plotting.base import GridspecPlot
    from models.plotting.kernel import KernelSubplot
except ModuleNotFoundError as e:
    raise e

cfg = Config(    
    mass_resolution=200,
    mass_max_value=1e12,
)

mg = DiscreteMassAxis(cfg)
mc = mg.bin_centers
ac = mg.particle_radii

R_coll = collision_cross_section(cfg, mg)

s = KernelSubplot(
    cfg, mg, R_coll,
    title=r"Collision Cross Section $\sigma^{coll}_{ij}$ [m$^2$]",
    # axis_scales=("log", "log", "lin"),
    cmap="Blues",
)

path_to_figures = Path(PATH_TO_FIGURES, "22")
os.makedirs(path_to_figures, exist_ok=True)
path_to_outfile = Path(path_to_figures, "collision_cross_section.pdf")

p = GridspecPlot([s])
p.render(
    save_plot=True,
    path_to_outfile=path_to_outfile
)
