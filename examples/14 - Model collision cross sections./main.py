import os
import sys
import matplotlib.pyplot as plt
try:
    sys.path.append(os.path.join("..", "..", "src"))
    from config import Config, PATH_TO_DARKMODE, PATH_TO_FIGURES
    from axis import DiscreteMassAxis
    from dust.collision_cross_section import collision_cross_section
except ModuleNotFoundError as e:
    raise e


# Define kernel configuration.
cfg = Config()

# Define discrete mass axis.
mg = DiscreteMassAxis(cfg)

# Setup pyplot figure.
FIGSIZE = (6, 5)
if cfg.mpl_dark_mode:
    plt.style.use(PATH_TO_DARKMODE)


if __name__ == "__main__":

    # Define collision cross section.
    sigma = collision_cross_section(cfg, mg)

    # Create plot.
    plt.pcolor(sigma)

    # Configure colorbar.
    plt.set_cmap("Reds")
    plt.colorbar()

    # Setup labels.
    plt.ylabel("$i$", rotation=0)
    plt.xlabel("$j$")
    plt.title(r"collision cross section $\sigma_{ij}$")

    # Save plot.
    os.makedirs("../../figures/21", exist_ok=True)
    path = os.path.join(PATH_TO_FIGURES, "21", "collision_cross_section.pdf")
    plt.savefig(path)
    plt.show()
    plt.close()
