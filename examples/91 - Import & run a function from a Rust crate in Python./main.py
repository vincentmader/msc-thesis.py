import numpy as np

try:
    import mc_coag
except ModuleNotFoundError:
    import sys
    sys.path.append("../../lib")
    try:
        import mc_coag
    except ModuleNotFoundError:
        print("Are you sure you have compiled `../../lib/mc_coag` ?")
        sys.exit()
from mc_coag.kernel import Kernel

# kernel = Kernel()
# print(kernel)
# K = kernel.entries()
# print(K)
# kernel = Kernel.from_config()
# print(kernel)

# ─────────────────────────────────────────────────────────────────────────────

# What do I have?
# - The complete kernel K_kij
# - The mass distribution N

# What do I want?
# - I want the kernel matrix K_kij to give to the integrator.
#   This has to be constructed in each step by sampling collisions.

m_min, m_max, N_m = 1e-10, 1e+10, 100
mass_grid = (m_min, m_max, N_m)

kernel = Kernel(mass_grid)
for _ in range(10):
    N = np.array([1, 0, 0, 1])  # TODO Use actual mass distribution.

    K = kernel.sample_from_mass_distribution(N)
    print(K)
