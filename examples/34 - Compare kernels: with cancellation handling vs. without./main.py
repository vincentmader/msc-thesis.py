import os, sys
import numpy as np
try:
    sys.path.append(os.path.join("..", "..", "src"))
    from axis import KernelAxis
    from config import Config
    from kernel import Kernel
    from visualization.base import GridspecPlot
    from visualization.kernel import KernelMassConservationSubplot, KernelSubplot
except ModuleNotFoundError as e:
    raise e

cfg_1 = Config(enable_cancellation_handling=False)
cfg_2 = Config(enable_cancellation_handling=True)
kernel_1, kernel_2 = Kernel(cfg_1), Kernel(cfg_2)
K_1, K_2 = kernel_1.K, kernel_2.K
K_diff = K_1 - K_2
K_equal = K_diff < 1e-14

mg = kernel_1.mg
mc = mg.bin_centers
ac = mg.particle_radii


def plot_1():
    p = GridspecPlot([
        KernelSubplot(
            mg, K_1,
            title="kernel $K_{kij}^{canc}$",
            scales=("log", "log", "lin"),
            symmetrized=True,
            cmap="bwr",
        ),
        KernelSubplot(
            mg, K_2,
            title="kernel $K_{kij}^{nocanc}$ with canc. handling",
            scales=("log", "log", "lin"),
            symmetrized=True,
            cmap="bwr",
        ),
        KernelSubplot(
            mg, np.abs(K_diff),
            title="abs($K_{kij}^{canc}-K_{kij}^{nocanc}$)",
            scales=("log", "log", "log"),
            symmetrized=True,
    )], add_slider=True)
    p.render()


def plot_2():
    p = GridspecPlot([
        KernelMassConservationSubplot(
            mg, K_1, axis=KernelAxis.Radius,
            title=r"kernel error $\Delta K^{canc}_{ij}=\sum_k m_k\cdot K_{kij}^{canc}$",
        ),
        KernelMassConservationSubplot(
            mg, K_2, axis=KernelAxis.Radius,
            title=r"kernel error $\Delta K^{nocanc}_{ij}=\sum_k m_k\cdot K_{kij}^{nocanc}$",
        ),
    ])
    p.render()


if __name__ == "__main__":
    plot_1()
    plot_2()
