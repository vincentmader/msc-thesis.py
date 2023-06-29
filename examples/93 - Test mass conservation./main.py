import os
import sys
import numpy as np
try:
    sys.path.append(os.path.join("..", "..", "src"))
    from axis import MassGrid
    from config import Config
    from kernel import Kernel
    from kernel import test_mass_conservation
    from visualization.kernel.v1.kernel_layer_plot import KernelLayerPlot
    from visualization.kernel.v1.interactive_kernel_layer_plot import InteractiveKernelLayerPlot
except ModuleNotFoundError as e:
    raise e


# Define configuration.
cfg = Config(
    mass_axis_scale="lin",
    mass_min_value=1,
    mass_max_value=51,
    mass_resolution=50,
    enable_coagulation=True,
    enable_fragmentation=False,
    enable_physical_cross_sections=False,
    enable_physical_relative_velocities=[],
    enable_fragmentation_variant=[
        # "naive/pulverization",
        # "mrn",
    ],
    enable_cancellation_handling=False,
)

# Setup pyplot figure.
FIGSIZE = (6, 5)
SLIDER_POSITION = [0.05, 0.3, 0.02, 0.4]

# Define discrete mass axis.
mg = MassGrid(cfg)
N_m = mg.N
mc = mg.grid_cell_centers
dm = mg.grid_cell_widths

# Define kernel.
kernel = Kernel(cfg)
K = kernel.K
K = np.array([0.5 * (K_k + K_k.T) for K_k in K])
# from kees_kernel import create_coag_kernel
# mgrain = mg.grid_cell_centers  # TODO Use mgrain instead of mg as well?
# K = create_coag_kernel(mgrain, R_coll)

# print(dms)
# for dm in dms:
#     print(dm)
# NOTE: The 2 print statements above do not print the exact same results!

# Count how many entries differ significantly from what they should be.
summa_2 = 0
a = []
for i in range(N_m):
    b = []
    for j in range(N_m):

        summa = 0
        for k in range(N_m - 1):  # TODO
            m_k = mc[k]
            # print(k, m_k)
            dm_k = dm[k]
            if cfg.mass_axis_scale == "lin" and np.abs(dm_k - 1) > 1e-14:
                raise Exception(
                    f"Grid spacing dm_k != 1 on linear grid, is this on purpose? (for {k=} -> {dm_k=})")

            summa += m_k * dm_k * K[k, i, j]

        machine_precision = abs(summa) < 1e-12
        b.append(machine_precision)
        # b.append(np.log(abs(summa)))
        # b.append(abs(summa))
        if machine_precision:
            summa_2 += 1
        print(f"{i=}\t{j=}\t{machine_precision}\t{summa=}")
    a.append(b)
print(f"\n{summa_2} = / {N_m**2} = {round(summa_2/N_m**2*100)} %")
b = np.array(a)

# Create plot & show it.
# plt.imshow(a, cmap="brg")
# plt.ylabel("$i$", rotation=0)
# plt.xlabel("$j$")

a = test_mass_conservation(kernel)

assert a.all() == b.all()

p = InteractiveKernelLayerPlot(
    [K], kernel_subplot_titles=["$K_{kij}$"],
    cmap_limits=(-1, 1),
)
p.show()

p = KernelLayerPlot(
    [a], kernel_subplot_titles=[
        r"$\sum_{k}\ K_{kij}\ m_k\ \Delta m_k\ \ (\overset{!}{=}0$ for mass conservation)"
    ],
)
p.show()
