import os, sys
try:
    sys.path.append(os.path.join("..", "..", "src"))
    from config import Config
    from kernel import Kernel
    from visualization.v1.interactive_kernel_layer_plot import InteractiveKernelLayerPlot
except ModuleNotFoundError as e:
    raise e


# Load configuration from `../../config.toml`.
cfg = Config(
    # mass_axis_scale="lin",
    # mass_resolution=50,
    # mass_min_value=2,
    # mass_max_value=52,

    # # mass_min_value=1e-4,
    # # mass_max_value=1e+4,

    # # enable_coagulation=True,
    # # enable_fragmentation=True,
    # enable_cancellation_handling=True,
    # enable_physical_collisions=False,
    # relative_velocity_components=[],

    # fragmentation_variant="mrn",
)

# Define kernel.
kernel = Kernel(cfg)
K = kernel.K

mg = kernel.mg


def main():

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
