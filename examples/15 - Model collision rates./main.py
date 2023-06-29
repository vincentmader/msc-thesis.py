import os
import sys

import matplotlib.pyplot as plt
try:
    sys.path.append(os.path.join("..", "..", "src"))
    from axis import MassGrid, RadialGrid
    from config import Config, PATH_TO_DARKMODE, PATH_TO_FIGURES
    from disk import Disk, DiskRegion
    from dust.collision_rate import collision_rate
except ModuleNotFoundError as e:
    raise e


# Define kernel configuration.
cfg = Config()

# Setup pyplot figure.
FIGSIZE = (6, 6)
if cfg.mpl_dark_mode:
    plt.style.use(PATH_TO_DARKMODE)

# Define discrete axis for radial distance from star, as well as for mass.
rg = RadialGrid(cfg)
mg = MassGrid(cfg)

# Define disk, the position of interest in it, & the disk properties there.
disk = Disk(cfg, rg, mg)
disk_region = DiskRegion(cfg, disk)


if __name__ == "__main__":

    # Define dust particle collision rate.
    R = collision_rate(cfg, disk, disk_region)

    # Create plot.
    plt.pcolor(R)
    plt.axis("scaled")

    # Define colorbar.
    plt.colorbar()
    plt.set_cmap("Reds")

    # Configure labels.
    plt.title("Collision rate $R_{coll}$")
    plt.ylabel("$i$", rotation=0)
    plt.xlabel("$j$")

    # Save plot.
    os.makedirs("../../figures/23", exist_ok=True)
    path = os.path.join(PATH_TO_FIGURES, "23", "collision_rate.pdf")
    plt.savefig(path)
    plt.show()
    plt.close()
