import os, sys
import numpy as np

sys.path.append(os.path.join("..", "..", "src"))
from config import Config
from config import PATH_TO_OUTFILES

from models.axis import DiscreteMassAxis
from models.kernel import Kernel
from models.kernel import SampledKernel

N_m = 100
rho_s = 0.5

N_half = int((N_m**2 + N_m) / 2)
# nr_of_samples = int(N_m * N_m / 2)
# nr_of_samples = int(rho_s * N_half) if rho_s != 1.0 else N_half
nr_of_samples = int(rho_s * N_half) if rho_s != 1.0 else N_m**2


cfg = Config(    
    enable_collision_sampling=False,
    mass_resolution=N_m,
)
K = Kernel(cfg)

cfg_2 = Config(
    enable_collision_sampling=True,
    mass_resolution=N_m,
    nr_of_samples=nr_of_samples,
)
mg = DiscreteMassAxis(cfg_2)
N = np.zeros(shape=(mg.N))
S = SampledKernel(cfg_2, N, W_ij=np.ones(shape=(N_m, N_m)))

from models.plotting.base import GridspecPlot
from models.plotting.evolution import EvolutionPlot, MassConservationPlot
from models.plotting.kernel import KernelSubplot, KernelMassConservationSubplot

GridspecPlot([
    KernelSubplot(
        cfg, mg, K.K_gain, 
        title=r"kernel gain $G_{kij}$, $K_{kij}\cdot R_{ij}$ [$kg m^3 s^{-1}$]",
        # axis_scales=(scale, scale, scale),
        # z_limits=z_limits,
    ),
    KernelSubplot(
        cfg, mg, -K.K_loss,
        title=r"kernel loss $L_{kij}$, $K_{kij}\cdot R_{ij}$ [$kg m^3 s^{-1}$]",
        # axis_scales=(scale, scale, scale),
        # z_limits=z_limits,
        # ylabel="",
        # symmetrized=True,
    ),
    KernelSubplot(
        cfg, mg, S.K_gain, 
        title=r"kernel gain $G_{kij}$, $K_{kij}\cdot R_{ij}$ [$kg m^3 s^{-1}$]",
        # axis_scales=(scale, scale, scale),
        # z_limits=z_limits,
        # symmetrized=symmetrized,
    ),
    KernelSubplot(
        cfg, mg, -S.K_loss,
        title=r"kernel loss $L_{kij}$, $K_{kij}\cdot R_{ij}$ [$kg m^3 s^{-1}$]",
        # axis_scales=(scale, scale, scale),
        # z_limits=z_limits,
        # ylabel="",
        # symmetrized=symmetrized,
        # axis_label_variant=axis_label_variant,
        # cmap=cmap,
    ),
], add_slider=True).render()

print("\nanother")
assert (K.K == S.K).all(), "Is N_sample = N_m**2 ?"
print("happy landing")
