import os, sys
try:
    sys.path.append(os.path.join("..", "..", "src"))
    from axis import DiscreteMassAxis
    from config import Config
    from kernel import Kernel
    from visualization.v1.interactive_kernel_layer_plot import InteractiveKernelLayerPlot
except ModuleNotFoundError as e:
    raise e

# Define kernel configuration.
cfg = Config(
    # Define mass axis.
    # On a linear grid, if we want to reach a mass grid spacing of exactly one,
    # we have to chose `mass_resolution = mass_max_value - mass_min_value`.
    mass_axis_scale="log",
    mass_min_value=1e-2,
    mass_max_value=1e+2,
    mass_resolution=50,
    # Define processes to include in the simulation.
    enable_coagulation=True,
    enable_fragmentation=False,
    enable_cancellation_handling=False,
    enable_physical_collisions=False,
    relative_velocity_components=[],
)

# Define discrete mass axis.
mg = DiscreteMassAxis(cfg)

# Define kernel.
kernel = Kernel(cfg)
K = kernel.K


def main():
    InteractiveKernelLayerPlot(
        [K],
        symmetrize_kernels=True,
        kernel_subplot_titles=["$K_{kij}$"],
        cmap_limits=(-1, 1),
    ).show()
