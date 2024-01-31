import os, sys
from pathlib import Path
from matplotlib.colors import LogNorm, Normalize
import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm
sys.path.append(os.path.join("..", "..", "src"))
from config import Config
from config import PATH_TO_FIGURES, PATH_TO_OUTFILES
from constants import SECONDS_PER_YEAR
from models.axis import DiscreteMassAxis
from models.solver import SolverV2
from functions.utils.dates import format_seconds_as_years

N_m = 50
SAMPLING_DENSITIES = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]

for rho_sample in SAMPLING_DENSITIES:
    nr_of_samples = int((N_m**2 + N_m) / 2 * rho_sample)

    cfg = Config(
        nr_of_samples=nr_of_samples
    )
