import os
import sys

import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.ticker import LogLocator, LogFormatterMathtext
try:
    sys.path.append(os.path.join("..", "..", "src"))
    from axis import DiscreteMassAxis, DiscreteRadialAxis
    from config import Config, PATH_TO_DARKMODE, PATH_TO_FIGURES
    from collision import collision_outcome_probabilities_from_maxwell_boltzmann 
    from collision import collision_outcome_probabilities_from_cutoff_velocity
    from disk import Disk, DiskRegion
    from dust.relative_velocity import relative_velocity
    from dust import particle_radius_from_mass
except ModuleNotFoundError as e:
    raise e


# Define kernel configuration.
cfg = Config()
rho_s = cfg.dust_particle_density

# Setup pyplot figure.
FIGSIZE = (6, 6)
if cfg.mpl_dark_mode:
    plt.style.use(PATH_TO_DARKMODE)

# Define discrete axis for radial distance from star, as well as for mass.
rg = DiscreteRadialAxis(cfg)
mg = DiscreteMassAxis(cfg)
mc = mg.grid_cell_centers
ac = particle_radius_from_mass(mc, rho_s)

# Define disk, the position of interest in it, & the disk properties there.
disk = Disk(cfg, rg, mg)
disk_region = DiskRegion(cfg, disk)

dv = relative_velocity(cfg, disk, disk_region)


def plot(title, P):

    # Create a grid of subplots (for pcolor plot & colorbar). 
    plt.figure()
    plt.set_cmap("Blues")
    gs = GridSpec(1, 2, width_ratios=[15, 1], wspace=0.)
    
    # Create the pcolor plot.
    ax1 = plt.subplot(gs[0])
    im = ax1.pcolor(ac, ac, P)
    plt.axis("scaled")
    
    # Set logarithmic scale for x & y axes.
    ax1.set_xscale("log")
    ax1.set_yscale("log")
    
    # Set logarithmic scale for x & y axis ticks.
    locator = LogLocator()
    ax1.xaxis.set_major_locator(locator)
    ax1.yaxis.set_major_locator(locator)
    
    # Set tick formatter with scientific notation.
    formatter = LogFormatterMathtext()
    ax1.xaxis.set_major_formatter(formatter)
    ax1.yaxis.set_major_formatter(formatter)
    
    # Add axis labels.
    ax1.set_ylabel("$a_i$ [m]", rotation=0)
    ax1.set_xlabel("$a_j$ [m]")
    
    # Create the colorbar.
    ax2 = plt.subplot(gs[1])
    plt.colorbar(im, cax=ax2)
    

if __name__ == "__main__":
    outcome_probabilities = {
        "cutoff_velocity": collision_outcome_probabilities_from_cutoff_velocity,
        "Maxwell-Boltzmann": collision_outcome_probabilities_from_maxwell_boltzmann,
    }

    for title, outcome_probability in outcome_probabilities.items():
        P_coag, P_frag = outcome_probability(cfg, dv)
        outcome_probabilities = {
            "P_coag": P_coag,
            "P_frag": P_frag,
        }

        for subtitle, P in outcome_probabilities.items():
            plot(title, P)

            # Show & save plot.
            filename = f"{subtitle} from {title}.pdf"
            os.makedirs("../../figures/16", exist_ok=True)
            path = os.path.join(PATH_TO_FIGURES, "16", filename)
            plt.savefig(path, bbox_inches='tight')
            plt.show()
            plt.close()
