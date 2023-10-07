import os, sys
import numpy as np
try:
    sys.path.append(os.path.join("..", "..", "src"))
    from axis import DiscreteMassAxis
    from config import Config
    from kernel import Kernel, SampledKernel
    from visualization.preset.p1 import plot_kernel_gain_loss, plot_kernel_error
except ModuleNotFoundError as e:
    raise e

N_m = 50
cfg = Config(
    mass_resolution=N_m,
    nr_of_samples=N_m**2,
)
mg = DiscreteMassAxis(cfg)

N = np.ones(shape=(N_m))

kernel_1 = Kernel(cfg)
kernel_2 = SampledKernel(cfg, N)

K = kernel_1.K
S = kernel_2.K

assert (K == S).all()

plot_kernel_gain_loss(cfg, mg, kernel_1, z_limits=(1e-20, 1), scale="log")
plot_kernel_gain_loss(cfg, mg, kernel_2, z_limits=(1e-20, 1), scale="log")
plot_kernel_error(cfg, mg, kernel_1, scale="log", z_limits=(1e-20, 1))
plot_kernel_error(cfg, mg, kernel_2, scale="log", z_limits=(1e-20, 1))
