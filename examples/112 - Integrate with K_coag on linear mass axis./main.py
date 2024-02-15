import os, sys
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
sys.path.append(os.path.join("..", "..", "src"))
from config import Config
from config import Config
from config import PATH_TO_FIGURES, PATH_TO_OUTFILES
from constants import AU
from functions.utils.dates import format_seconds_as_years
from models.axis import DiscreteTimeAxis
from models.disk import Disk, DiskRegion
from models.kernel import Kernel
from models.solver import SolverV2

path_to_outfiles = Path(PATH_TO_OUTFILES, "data", "112")
path_to_figures  = Path(PATH_TO_FIGURES, "112")
os.makedirs(path_to_figures, exist_ok=True)

dt  = 1
N_m = 100
N_t = 200
cfg = Config(
    mass_min_value = 1,
    mass_max_value = 1e6,
    mass_resolution = N_m,
    enable_coagulation = True,
    enable_fragmentation = False,
    enable_collision_sampling = False, # note: `True` can not be use in the context of this module file.
)

kernel = Kernel(cfg, R_coag=np.ones([N_m]*2))
K = kernel.K
mg = kernel.mg
mc = mg.bin_centers
tg = DiscreteTimeAxis(cfg)
tc = tg.bin_centers

N_dust = np.array([0.0] * N_m)
N_dust[0] = 1.0
N_dust_vs_t = []

for t in range(N_t):

    N_dust_vs_t.append(N_dust.copy())

    N_i = N_dust[None, :, None]
    N_j = N_dust[None, None, :]
    dNdt = (K[:, :, :] * N_i * N_j).sum(axis=2).sum(axis=1)
    N_dust += dNdt * dt

RELEVANT_TIMES = np.logspace(np.log10(20), np.log10(N_t-1), 10)
print(RELEVANT_TIMES)
# RELEVANT_TIMES = np.array([0, 0.4, 0.2, 0.6, 0.8]) * N_m
# RELEVANT_TIMES = [0] + [int(i) for i in RELEVANT_TIMES] + [N_m-1]
RELEVANT_TIMES = [int(i) for i in RELEVANT_TIMES]
colors = plt.cm.Reds(np.linspace(0.2, 1, len(RELEVANT_TIMES)))

plt.figure(figsize=(10, 4))

for p, t in enumerate(RELEVANT_TIMES):
    t_str = format_seconds_as_years(tc[t])

    N_dust = N_dust_vs_t[t]
    plt.loglog(mc, mc * N_dust, color=colors[p], label=f"t={t_str}")

plt.grid(True)
plt.ylim(1e-5, 10)
# plt.gca().set_ylim(bottom=1e-10)
plt.legend(loc="best")
plt.xlabel("Dust Particle Mass $m_i^c$ [kg]")
plt.ylabel(r"Dust Density $\rho_i^d$ [kg m$^{-3}$]")
plt.savefig(Path(path_to_figures, "euler_explicit_integration_on_linear_mass_axis.pdf"))

plt.show()
plt.close()
