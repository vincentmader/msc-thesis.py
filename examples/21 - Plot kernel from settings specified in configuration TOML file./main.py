from pprint import pprint
import os
import sys
try:
    sys.path.append(os.path.join("..", "..", "src"))
    from config import Config
    from kernel import Kernel
    from visualization.kernel.v1.interactive_kernel_layer_plot import InteractiveKernelLayerPlot
    from visualization.kernel.v2.gain_vs_loss import KernelGainVsLossPlot
    from visualization.kernel.v2.mass_conservation import KernelMassConservationPlot
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
    # enable_physical_cross_sections=False,
    # enable_physical_relative_velocities=[],

    # fragmentation_variant="mrn",
)
pprint(cfg.__dict__)

# Define kernel.
kernel = Kernel(cfg)
K = kernel.K

mg = kernel.mg


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

    # Plot kernel gain & loss terms separately.
    p = KernelGainVsLossPlot(kernel)
    p.show()

    # Plot `\sum_{ij} m_k \Delta m_k K_kij`.
    p = KernelMassConservationPlot(cfg, mg, K)
    p.show()
