import os, sys
from pathlib import Path
try:
    sys.path.append(os.path.join("..", "..", "src"))
    from models.axis import DiscreteMassAxis, DiscreteRadialAxis, AxisLabelVariant
    from collision import collision_rate
    from config import Config, PATH_TO_FIGURES
    from disk import Disk, DiskRegion
    from visualization.base import GridspecPlot
    from visualization.kernel import KernelSubplot
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
disk = Disk(cfg, rg, mg)
disk_region = DiskRegion(cfg, disk)
R_coll = collision_rate(cfg, disk, disk_region)

path_to_figures = Path(PATH_TO_FIGURES, "23")
os.makedirs(path_to_figures, exist_ok=True)
path_to_outfile = Path(path_to_figures, "collision_rate.pdf")

p = GridspecPlot([
    KernelSubplot(
        cfg, mg, R_coll,
        title="collision rate $R_{coll}$",
        # ^ TODO Add units.
        axis_label_variant=AxisLabelVariant.Radius,
        # axis_scales=("log", "log", "lin"),
        cmap="Blues",
    )
]).render(
    save_plot=True,
    path_to_outfile=path_to_outfile
)
