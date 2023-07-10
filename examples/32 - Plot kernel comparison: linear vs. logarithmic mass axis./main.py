import os
import sys
import numpy as np
try:
    sys.path.append(os.path.join("..", "..", "src"))
    from axis import DiscreteMassAxis
    from config import Config
    from kernel import Kernel
    from visualization.kernel.v1.interactive_kernel_layer_plot import InteractiveKernelLayerPlot
except ModuleNotFoundError as e:
    raise e


# Define kernel configuration.
cfg = Config(
    enable_coagulation=True,
    enable_fragmentation=False,
    enable_physical_collisions=False,
    relative_velocity_components=[],
)


if __name__ == "__main__":

    # Define list of kernel subplot titles.
    kernel_subplot_titles = [
        "$K_{kij}$",
        r"$\frac{1}{2}(K_{kij} + K_{kij}^T)$",
    ]

    for scale in ["lin", "log"]:
        cfg.mass_axis_scale = scale

        # Define discrete mass axis.
        mg = DiscreteMassAxis(cfg)

        # Define kernel.
        kernel = Kernel(cfg)
        K_1 = kernel.K
        K_2 = np.array([0.5 * (K_k + K_k.T) for K_k in K_1])

        # Define list of kernels to plot.
        kernels = [K_1, K_2]

        # Create plot & show it.
        p = InteractiveKernelLayerPlot(
            kernels,
            kernel_subplot_titles=kernel_subplot_titles,
            cmap_limits=(-1, 1),
        )
        p.show()
