import os, sys
from pathlib import Path
try:
    sys.path.append(os.path.join("..", "..", "src"))
    from config import Config, PATH_TO_FIGURES
    from functions.dust.collision import collision_rate
    from functions.dust.collision import collision_outcome_probabilities_from_cutoff_velocity
    from functions.dust.collision import collision_outcome_probabilities_from_maxwell_boltzmann
    from functions.dust.relative_velocity import relative_velocity
    from models.axis import DiscreteMassAxis, DiscreteRadialAxis
    from models.disk import Disk, DiskRegion
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

disk = Disk(cfg, rg, mg)
disk_region = DiskRegion(cfg, disk)
dv = relative_velocity(cfg, disk, disk_region)

R_coll = collision_rate(cfg, disk, disk_region)

fs = [
    collision_outcome_probabilities_from_cutoff_velocity, 
    collision_outcome_probabilities_from_maxwell_boltzmann,
]
for label_f, f in zip(["cutoff", "MB"], fs):
    P_coag, P_frag = f(cfg, dv)
    R_coag = R_coll * P_coag
    R_frag = R_coll * P_frag

    for label_R, R in zip(["coag", "frag"], [R_coag, R_frag]):
    
        s = PcolorMatrixSubplot(
            ac, ac, R,
            title=label_R + ". rate $R_{" + label_R + "}$",
            # ^ TODO Add units.
            xlabel="particle radius $a_j$ [m]",
            ylabel="particle radius $a_i$ [m]",
            z_limits=(1e-30, 1e5), # TODO Determine dynamically?
            cmap="Blues",
            # axis_scales=("log", "log", "log"), # TODO Define dynamically.
        )
        
        path_to_figures = Path(PATH_TO_FIGURES, "25")
        os.makedirs(path_to_figures, exist_ok=True)
        path_to_outfile = Path(path_to_figures, f"R_{label_R}_from_{label_f}.pdf")
        
        GridspecPlot([s]).render(
            save_plot=True,
            path_to_outfile=path_to_outfile
        )
