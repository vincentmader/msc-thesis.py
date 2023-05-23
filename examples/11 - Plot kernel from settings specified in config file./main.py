import os
import sys
try:
    sys.path.append(os.path.join("..", "..", "src"))
    from config import Config
    from disk import MassGrid, Disk, RadialGrid, DiskRegion
    from kernel import Kernel
    from kernel.collision_rate import collision_rate
    from visualization.kernel.interactive_kernel_layer_plot import InteractiveKernelLayerPlot
except ModuleNotFoundError as e:
    raise e


# Load configuration from `../../config.toml`.
cfg = Config()

# Define discrete axis for radial distance from star, as well as for mass.
rg = RadialGrid(cfg)
mg = MassGrid(cfg)

# Define disk, the position of interest in it, & the disk properties there.
disk = Disk(cfg, rg, mg)
disk_region = DiskRegion(cfg, disk)
R_coll = collision_rate(cfg, disk, disk_region)

# Define kernel.
kernel = Kernel(cfg)
K = kernel.K(mg, R_coll)


if __name__ == "__main__":

    # Define list of kernels to plot.
    kernels = [K]

    # Create plot & show it.
    p = InteractiveKernelLayerPlot(
        kernels,
        symmetrize_kernels=True,
        kernel_subplot_titles=["$K_{kij}$"],
    )
    p.show()
