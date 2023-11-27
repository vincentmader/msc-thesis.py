import os, sys
import numpy as np
try:
    sys.path.append(os.path.join("..", "..", "src"))
    from config import Config
    from models.axis import DiscreteMassAxis
    from models.kernel import Kernel, SampledKernel
    from visualization.preset.p1 import plot_kernel_gain_loss, plot_kernel_error
except ModuleNotFoundError as e:
    raise e

# Define configuration.
N_m = 50
cfg = Config(
    mass_resolution=N_m,
    nr_of_samples=N_m**2,
)
mg = DiscreteMassAxis(cfg)

# Define number of particles per bin. (-> Here: 1)
N = np.ones(shape=(N_m))

# Define kernels, both with & without sampling.
kernel_1 = Kernel(cfg)
K = kernel_1.K
kernel_2 = SampledKernel(cfg, N)
S = kernel_2.K

# Assert equality.
assert (K == S).all()
# NOTE: This assertion actually works perfectly fine!
#    -> The sampled kernel does exactly equal the unsampled one,
#       if we let `N_samples = N_m^2`.
# NOTE: This does not mean that everything works though!
#       The assertion will fail in the solver.
#    -> Why? Because `N` will have entries equal to zero!
#    -> Then `P_ij = 0` as well.
#    -> This means that some collisions will NEVER be sampled,
#       regardless of how many points I select at random.
# NOTE: The solution to this should be quite simple:
#    -> Set `P_ij = 1e-99` instead of 0 (or some other tiny value).
#    -> This should affect the result ONLY in the limit of `N_samples = N_m^2`.

# Plot the kernels, & their behavior re: mass error.
plot_kernel_gain_loss(cfg, mg, kernel_1, z_limits=(1e-20, 1), scale="log")
plot_kernel_gain_loss(cfg, mg, kernel_2, z_limits=(1e-20, 1), scale="log")
plot_kernel_error(cfg, mg, kernel_1, scale="log", z_limits=(1e-20, 1))
plot_kernel_error(cfg, mg, kernel_2, scale="log", z_limits=(1e-20, 1))
