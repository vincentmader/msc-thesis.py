import os
import sys
try:
    sys.path.append(os.path.join("..", "..", "src"))
    from config import Config
    from disk import MassGrid
    from kernel import Kernel
    from visualization.kernel.v1.interactive_kernel_layer_plot import InteractiveKernelLayerPlot
except ModuleNotFoundError as e:
    raise e


# Define kernel configuration.
cfg = Config(
    # Define mass axis.
    # On a linear grid, if we want to reach a mass grid spacing of exactly one,
    # we have to chose `mass_resolution = mass_max_value - mass_min_value`.
    mass_axis_scale="lin",
    mass_min_value=1,
    mass_max_value=51,
    mass_resolution=50,
    # Define processes to include in the simulation.
    enable_coagulation=False,
    enable_fragmentation=True,
    enable_cancellation_handling=False,
    enable_physical_cross_sections=False,
    enable_physical_relative_velocities=[],
    enable_fragmentation_variant=[
        # "naive/pulverization",
        "mrn",
    ],
)

# Define discrete mass axis.
mg = MassGrid(cfg)

# Define kernel.
kernel = Kernel(cfg)
K = kernel.K


if __name__ == "__main__":

    # Define list of kernels to plot.
    kernels = [K]

    # Create plot & show it.
    p = InteractiveKernelLayerPlot(
        kernels,
        symmetrize_kernels=True,
        kernel_subplot_titles=["$K_{kij}$"],
        cmap_limits=(-1, 1),
    )
    p.show()
