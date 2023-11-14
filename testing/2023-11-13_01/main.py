import numpy as np

PATH_TO_FILES = "/Users/vinc/Desktop/exports"
# ^ Redefine this to reflect the location of the export files on your machine.

# Load mass grid (values at bin center).
path = f"{PATH_TO_FILES}/m_k.txt"
m_k = np.loadtxt(path)

# Load bin widths.
path = f"{PATH_TO_FILES}/dm_k.txt"
dm_k = np.loadtxt(path)

# Load particle mass distribution values. 
# - These are the values at `i_t = 100` -> `t ~ 51y`.
# - Note: `N_k = n_k * dm_k` and `rho = n_k * m_k * dm_k`
path = f"{PATH_TO_FILES}/n_k.txt"
n_k = np.loadtxt(path)

# Define number of bins in mass grid:
N_m = m_k.shape[0]

# Load & reshape complete kernel.
path = f"{PATH_TO_FILES}/K_kij_complete.txt"
K_kij_complete = np.loadtxt(path)
K_kij_complete = K_kij_complete.reshape((N_m, N_m, N_m))

# Load & reshape sampled kernel.
path = f"{PATH_TO_FILES}/K_kij_sampled.txt"
K_kij_sampled = np.loadtxt(path)
K_kij_sampled = K_kij_sampled.reshape((N_m, N_m, N_m))
