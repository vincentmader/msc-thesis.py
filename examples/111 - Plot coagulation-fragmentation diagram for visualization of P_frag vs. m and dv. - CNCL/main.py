import os, sys
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
sys.path.append(os.path.join("..", "..", "src"))
from config import Config
from config import PATH_TO_FIGURES, PATH_TO_OUTFILES
from models.solver import SolverV2
from functions.dust.collision import collision_outcome_probabilities_from_maxwell_boltzmann 
from functions.dust.collision import collision_outcome_probabilities_from_cutoff_velocity
from functions.dust.relative_velocity import relative_velocity
from models.axis import DiscreteMassAxis, DiscreteRadialAxis
from models.disk import Disk, DiskRegion
from models.plotting.base import GridspecPlot, PcolorMatrixSubplot

path_to_outfiles = Path(PATH_TO_OUTFILES, "data", "104")
path_to_figures  = Path(PATH_TO_FIGURES, "104")

# Define kernel configuration.
cfg = Config(
    mass_resolution=200,
    mass_max_value=1e12,
)

# Define discrete axis for radial distance from star, as well as for mass.
rg = DiscreteRadialAxis(cfg)
mg = DiscreteMassAxis(cfg)
mc = mg.bin_centers
ac = mg.particle_radii

# Define disk, the position of interest in it, & the disk properties there.
disk = Disk(cfg, rg, mg)
disk_region = DiskRegion(cfg, disk)

dv = relative_velocity(cfg, disk, disk_region)

outcome_probabilities = {
    "cutoff_velocity": collision_outcome_probabilities_from_cutoff_velocity,
    "Maxwell-Boltzmann": collision_outcome_probabilities_from_maxwell_boltzmann,
}

for label, outcome_probability in outcome_probabilities.items():

    P_coag, P_frag = outcome_probability(cfg, dv)
    outcome_probabilities = {
        "P_coag": P_coag,
        "P_frag": P_frag,
    }


v_frag = cfg.fragmentation_velocity
N_m = 100
for i in range(N_m):
    for j in range(N_m):
        P_f = (3/2 * (v_frag/dv[i, j])**2 + 1) * np.exp(-3/2 * (v_frag/dv[i, j])**2) # See "2022 Stammler & Birnstiel"
        P_c = 1 - P_f # NOTE: Bouncing is neglected here.
