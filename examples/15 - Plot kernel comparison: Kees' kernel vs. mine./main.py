import os
import sys
import numpy as np
try:
    sys.path.append(os.path.join("..", "..", "src"))
    from config import Config
    from disk import MassGrid
    from kees_kernel import create_coag_kernel
    from kernel import Kernel
    from visualization.kernel.interactive_kernel_layer_plot import InteractiveKernelLayerPlot
except ModuleNotFoundError as e:
    raise e


# Define kernel configuration.
cfg = Config(
    mass_axis_scale="lin",
    mass_resolution=50,
    mass_min_value=2,
    mass_max_value=52,
    enable_coagulation=True,
    enable_fragmentation=False,
    enable_cancellation_handling=False,
    enable_physical_cross_sections=False,
    enable_physical_relative_velocities=[],
    # mass_open_boundary=True,
)

# Define discrete mass axis.
mg = MassGrid(cfg)
N_m = mg.N

# Define collision rate.
# Here: Set to 1 for simplicity.
R_coll = np.ones(shape=[N_m]*2)

# Define kernel according to definition in `../../src/kernel/__init__.py`.
kernel = Kernel(cfg)
Kkij_vinc = kernel.K
Kkij_vinc = np.array([0.5 * (K_k + K_k.T) for K_k in Kkij_vinc])

# Define kernel according to Kees' definition in `./kees_kernel.py`.
mgrain = mg.grid_cell_boundaries()[:-1]  # TODO Use mgrain instead of mg as well?
Kkij_kees = create_coag_kernel(mgrain, R_coll)  # NOTE: R_coll was originally named Cij

# Define comparison metrics.
Kkij_equal = (Kkij_vinc - Kkij_kees) < 1e-14
    # Kkij_diff = Kkij_kees - Kkij_vinc
    # Kkij_v2k = (Kkij_vinc - Kkij_kees) / Kkij_kees * 100
    # Kkij_k2v = (Kkij_kees - Kkij_vinc) / Kkij_vinc * 100
    # Kkij_log_v2k = np.log(Kkij_v2k)
    # Kkij_log_k2v = np.log(Kkij_k2v)


if __name__ == "__main__":

    # Define list of kernels to plot.
    kernels = [
        Kkij_kees,
        Kkij_vinc,
        Kkij_equal,
        # Kkij_diff,
        # Kkij_k2v,
        # Kkij_v2k,
    ]

    # Define list of kernel subplot titles.
    kernel_subplot_titles = [
        "$K_{kees}$",
        "$K_{vinc}$",
        "$K_{kees}$ == $K_{vinc}$",
    ]

    # Create plot & show it.
    p = InteractiveKernelLayerPlot(
        kernels,
        kernel_subplot_titles=kernel_subplot_titles,
        cmap_limits=(-1, 1),
    )
    p.show()

    # ─────────────────────────────────────────────────────────────────────────

    # Calculate number of kernel entries differing sigificantly from Kees'.
    count_of_different_entries = 0
    count_of_significantly_different_entries = 0
    for k in range(N_m):
        Kv = Kkij_vinc[k]
        Kv = 0.5 * (Kv + Kv.T)
        Kk = Kkij_kees[k]
        for i in range(N_m):
            for j in range(N_m):
                if Kv[i,j] == Kk[i,j]:
                    continue
                count_of_different_entries += 1
                rel_diff = (Kk[i,j]-Kv[i,j])/Kv[i,j]
                if np.abs(rel_diff) > 1e-12:
                    # print(k, i, j, "\tKk=", Kk[i,j], "\tKv=", Kv[i,j], "\t", rel_diff * 100, "%")
                    count_of_significantly_different_entries += 1
    print(f"{count_of_different_entries=}")
    print(f"{count_of_significantly_different_entries=}")
