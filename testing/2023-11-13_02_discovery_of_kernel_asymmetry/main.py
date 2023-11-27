#!/usr/bin/env python3
import sys
sys.path.append("../../src")
import numpy as np
from config import Config
from models.axis import DiscreteMassAxis
from models.kernel import Kernel

# ─────────────────────────────────────────────────────────────────────────────

cfg = Config(
    enable_coagulation=True,
    enable_fragmentation=True,
)
mg = DiscreteMassAxis(cfg)

# ─────────────────────────────────────────────────────────────────────────────

N_included = 20
ijs = [(i, j) for i in range(mg.N) for j in range(20)]
# NOTE: I am deliberately NOT including all `(j,i)` for each `(i,j)`
#    -> Set `N_included = mg.N`, and the error disappears.
# NOTE: If this doesn't work: Asymmetry is handled in the kernel definition.

# ─────────────────────────────────────────────────────────────────────────────

R = np.ones(shape=(mg.N, mg.N))
kernel = Kernel(
    cfg, 
    ijs=ijs, 
    # R_coag=R, R_frag=R,
)
R = kernel.R_coag + kernel.R_frag

# ─────────────────────────────────────────────────────────────────────────────

scale = "log"
z_limits = (1e-20, 1)

from functions.plotting.preset.p1 import plot_kernel_gain_loss
plot_kernel_gain_loss(cfg, mg, kernel, scale, z_limits)

from functions.plotting.preset.p1 import plot_kernel_error
plot_kernel_error(cfg, mg, kernel, R, scale, z_limits, symmetrized=False)
