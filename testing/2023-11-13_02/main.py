#!/usr/bin/env python3
import sys
sys.path.append("../../src")

from axis import DiscreteMassAxis
from config import Config
from kernel import Kernel

scale = "log"
z_limits = (1e-20, 1)

cfg = Config(
    # mass_resolution=6,
    # mass_min_value=1e-8,
    # mass_max_value=1e-7,
)
mg = DiscreteMassAxis(cfg)

ijs = [(i, j) for i in range(mg.N) for j in range(20)]

kernel = Kernel(cfg, ijs=ijs)
R = kernel.R_coag + kernel.R_frag

# from visualization.preset.p1 import plot_kernel_gain_loss
# plot_kernel_gain_loss(cfg, mg, kernel, scale, z_limits)

from visualization.preset.p1 import plot_kernel_error
plot_kernel_error(cfg, mg, kernel, R, scale, z_limits, symmetrized=False)
