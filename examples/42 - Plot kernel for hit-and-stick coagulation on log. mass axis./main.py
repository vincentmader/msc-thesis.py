import os, sys
try:
    sys.path.append(os.path.join("..", "..", "src"))
    from axis import KernelAxis
    from config import Config
    from kernel import Kernel
    from visualization.base import GridspecPlot
    from visualization.kernel import KernelSubplot
except ModuleNotFoundError as e:
    raise e

cfg = Config(
    enable_fragmentation=False,
    enable_physical_collisions=False,
    relative_velocity_components=[],
)
kernel = Kernel(cfg)
mg = kernel.mg


def plot_1():
    GridspecPlot([
        KernelSubplot(
            mg, kernel.K,
            title="kernel $K_{kij}$",
            scales=("lin", "lin", "lin"),
            symmetrized=True,
            axis=KernelAxis.Bin,
            cmap="bwr",
            z_limits=(-1, 1)
        )
    ], add_slider=True).render()


def plot_2():
    GridspecPlot([
        KernelSubplot(
            mg, kernel.K_gain, 
            title="kernel gain contribution $G_{kij}$",
            scales=("lin", "lin", "lin"),
            symmetrized=True,
            axis=KernelAxis.Bin,
        ),
        KernelSubplot(
            mg, -kernel.K_loss,
            title="kernel loss contribution $L_{kij}$",
            scales=("lin", "lin", "lin"),
            symmetrized=True,
            axis=KernelAxis.Bin,
        ),
    ], add_slider=True).render()

if __name__ == "__main__":
    plot_1()
    plot_2()
