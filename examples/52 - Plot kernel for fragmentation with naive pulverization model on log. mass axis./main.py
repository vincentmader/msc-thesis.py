import os, sys
try:
    sys.path.append(os.path.join("..", "..", "src"))
    from axis import KernelAxis
    from config import Config
    from kernel import Kernel
    from visualization.base import GridspecPlot
    from visualization.kernel.kernel import KernelSubplot
except ModuleNotFoundError as e:
    raise e

cfg = Config(
    enable_coagulation=False,
    enable_fragmentation=True,
    enable_cancellation_handling=True,
    enable_physical_collisions=False,
    relative_velocity_components=[],
    fragmentation_variant="naive/pulverization",
)
kernel = Kernel(cfg)
K, mg = kernel.K, kernel.mg


def plot_1():
    GridspecPlot([
        KernelSubplot(
            mg, kernel.K,
            title="kernel $K_{kij}$",
            scales=("log", "log", "lin"),
            symmetrized=True,
            axis=KernelAxis.Radius,
            cmap="bwr",
            z_limits=(-1, 1)
        )
    ], add_slider=True).render()


def plot_2():
    GridspecPlot([
        KernelSubplot(
            mg, kernel.K_gain, 
            title="kernel gain contribution $G_{kij}$",
            scales=("log", "log", "log"),
            symmetrized=True,
            axis=KernelAxis.Radius,
        ),
        KernelSubplot(
            mg, -kernel.K_loss,
            title="kernel loss contribution $L_{kij}$",
            scales=("log", "log", "log"),
            symmetrized=True,
            axis=KernelAxis.Radius,
        ),
    ], add_slider=True).render()


if __name__ == "__main__":
    plot_1()
    plot_2()
