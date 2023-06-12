import os
import sys
try:
    sys.path.append(os.path.join("..", "..", "src"))
    from config import Config
    from kernel import Kernel
    from visualization.kernel.v1.interactive_kernel_layer_plot import InteractiveKernelLayerPlot
    from visualization.kernel.v2.kernel_gain_vs_loss import KernelGainVsLossPlot
    from visualization.kernel.v2.mass_conservation import KernelMassConservationPlot
except ModuleNotFoundError as e:
    raise e


# Load configuration from `../../config.toml`.
cfg = Config()

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

    # Plot kernel gain & loss terms separately.
    p = KernelGainVsLossPlot(kernel)
    p.show()

    # Plot `\sum_{ij} m_k \Delta m_k K_kij`.
    p = KernelMassConservationPlot(kernel)
    p.show()
