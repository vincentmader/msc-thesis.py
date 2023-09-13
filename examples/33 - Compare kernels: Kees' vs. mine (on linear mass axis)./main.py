import os
import sys

import numpy as np
try:
    sys.path.append(os.path.join("..", "..", "src"))
    from config import Config
    from dust import particle_radius_from_mass
    from kees_kernel import create_coag_kernel
    from kernel import Kernel
    from visualization.kernel.pcolor_matrix_subplot import PcolorMatrixSubplot
    from visualization.kernel.gridspec_plot import GridspecPlot
except ModuleNotFoundError as e:
    raise e


cfg = Config(
    mass_axis_scale="lin",
    mass_resolution=50,
    mass_min_value=2,
    mass_max_value=52,

    # mass_min_value=1e-4,
    # mass_max_value=1e+4,

    enable_coagulation=True,
    enable_fragmentation=False,
    enable_cancellation_handling=False,
    enable_physical_collisions=False,
    relative_velocity_components=[],
)

kernel_vinc = Kernel(cfg)
K_vinc = kernel_vinc.K

mg = kernel_vinc.mg
mc = mg.grid_cell_centers
rho_s = cfg.dust_particle_density
ac = particle_radius_from_mass(mc, rho_s)

N_m = mg.N
R_coll = np.ones(shape=[N_m] * 2)

K_kees = create_coag_kernel(mc, R_coll)

i = np.linspace(0, N_m, N_m)

K_diff = K_vinc - K_kees
K_equal = K_diff < 1e-14
# Kkij_v2k = (Kkij_vinc - Kkij_kees) / Kkij_kees * 100
# Kkij_k2v = (Kkij_kees - Kkij_vinc) / Kkij_vinc * 100
# Kkij_log_v2k = np.log(Kkij_v2k)
# Kkij_log_k2v = np.log(Kkij_k2v)

s1 = PcolorMatrixSubplot(
    # ac, ac, K_vinc,
    i, i, K_vinc,
    title="$K_{kij}^{vinc}$",
    xlabel="bin index $j$",
    ylabel="bin index $i$",
    scales=("lin", "lin", "lin"),
    symmetrized=True,
)
s2 = PcolorMatrixSubplot(
    # ac, ac, K_kees,
    i, i, K_kees,
    title="$K_{kij}^{kees}$",
    xlabel="bin index $j$",
    scales=("lin", "lin", "lin"),
    symmetrized=True,
)
s3 = PcolorMatrixSubplot(
    # ac, ac, K_kees,
    i, i, K_equal,
    title="$K_{kij}^{vinc}=K_{kij}^{kees}$",
    xlabel="bin index $j$",
    scales=("lin", "lin", "lin"),
    symmetrized=True,
)
p = GridspecPlot([s1, s2, s3], add_slider=True)
p.render()

# ═════════════════════════════════════════════════════════════════════════════

# count_of_different_entries = 0
# count_of_significantly_different_entries = 0
# for k in range(N_m):
#     Kv = Kkij_vinc[k]
#     Kv = 0.5 * (Kv + Kv.T)
#     Kk = Kkij_kees[k]
#     for i in range(N_m):
#         for j in range(N_m):
#             if Kv[i, j] == Kk[i, j]:
#                 continue
#             count_of_different_entries += 1
#             rel_diff = (Kk[i, j] - Kv[i, j]) / Kv[i, j]
#             if np.abs(rel_diff) > 1e-12:
#                 # print(k, i, j, "\tKk=", Kk[i,j], "\tKv=", Kv[i,j], "\t", rel_diff * 100, "%")
#                 count_of_significantly_different_entries += 1
# print(f"{count_of_different_entries=}")
# print(f"{count_of_significantly_different_entries=}")

# # Plot `\sum_{ij} m_k \Delta m_k K_kij`.
# print("Showing mass conservation test for `K_kij_vinc`...")
# p = KernelMassConservationPlot(cfg, mg, Kkij_vinc)
# p.show()

# # Plot `\sum_{ij} m_k \Delta m_k K_kij`.
# print("Showing mass conservation test for `K_kij_kees`...")
# p = KernelMassConservationPlot(cfg, mg, Kkij_kees)
# p.show()

from visualization.v1.mass_conservation import KernelMassConservationPlot

# Plot `\sum_{ij} m_k \Delta m_k K_kij`.
p = KernelMassConservationPlot(cfg, mg, K_vinc)
p.show()

# Plot `\sum_{ij} m_k \Delta m_k K_kij`.
p = KernelMassConservationPlot(cfg, mg, K_kees)
p.show()
