import os, sys
from pathlib import Path
import matplotlib.pyplot as plt
try:
    sys.path.append(os.path.join("..", "..", "src"))
    from axis import DiscreteMassAxis, DiscreteRadialAxis
    from collision import collision_outcome_probabilities_from_maxwell_boltzmann 
    from collision import collision_outcome_probabilities_from_cutoff_velocity
    from config import Config, PATH_TO_FIGURES
    from disk import Disk, DiskRegion
    from dust.relative_velocity import relative_velocity
    from visualization.base import GridspecPlot, PcolorMatrixSubplot
except ModuleNotFoundError as e:
    raise e

# Define kernel configuration.
cfg = Config()

# Define discrete axis for radial distance from star, as well as for mass.
rg = DiscreteRadialAxis(cfg)
mg = DiscreteMassAxis(cfg)
mc = mg.grid_cell_centers
ac = mg.particle_radii

# Define disk, the position of interest in it, & the disk properties there.
disk = Disk(cfg, rg, mg)
disk_region = DiskRegion(cfg, disk)

dv = relative_velocity(cfg, disk, disk_region)

outcome_probabilities = {
    "cutoff_velocity": collision_outcome_probabilities_from_cutoff_velocity,
    "Maxwell-Boltzmann": collision_outcome_probabilities_from_maxwell_boltzmann,
}


def plot_together(title, P_coag, P_frag):
    s1 = PcolorMatrixSubplot(
        ac, ac, P_coag,
        title="dust particle coagulation probability $P_{coag}$",
        xlabel="particle radius $a_j$ [m]",
        ylabel="particle radius $a_i$ [m]",
        scales=("log", "log", "lin"),
        cmap="Blues",
    )
    s2 = PcolorMatrixSubplot(
        ac, ac, P_frag,
        title="dust particle fragmentation probability $P_{frag}$",
        xlabel="particle radius $a_j$ [m]",
        ylabel="particle radius $a_i$ [m]",
        scales=("log", "log", "lin"),
        cmap="Blues",
    )
    
    path_to_figures = Path(PATH_TO_FIGURES, "24")
    os.makedirs(path_to_figures, exist_ok=True)
    filename = f"collision_outcome_probabilities_from_{title}.pdf"
    path_to_outfile = Path(path_to_figures, filename)
    
    p = GridspecPlot([s1, s2])
    p.render(
        save_plot=True,
        path_to_outfile=path_to_outfile
    )

def plot_separately(title, P_coag, P_frag):
    s = PcolorMatrixSubplot(
        ac, ac, P_coag,
        title="dust particle coagulation probability $P_{coag}$",
        xlabel="particle radius $a_j$ [m]",
        ylabel="particle radius $a_i$ [m]",
        scales=("log", "log", "lin"),
        cmap="Blues",
    ) 
    
    path_to_figures = Path(PATH_TO_FIGURES, "24")
    os.makedirs(path_to_figures, exist_ok=True)
    filename = f"coagulation_probability_from_{title}.pdf"
    path_to_outfile = Path(path_to_figures, filename)
    
    p = GridspecPlot([s])
    p.render(
        save_plot=True,
        path_to_outfile=path_to_outfile
    )

    s = PcolorMatrixSubplot(
        ac, ac, P_frag,
        title="dust particle fragmentation probability $P_{frag}$",
        xlabel="particle radius $a_j$ [m]",
        ylabel="particle radius $a_i$ [m]",
        scales=("log", "log", "lin"),
        cmap="Blues",
    )

    path_to_figures = Path(PATH_TO_FIGURES, "24")
    os.makedirs(path_to_figures, exist_ok=True)
    filename = f"fragmentation_probability_from_{title}.pdf"
    path_to_outfile = Path(path_to_figures, filename)
    
    p = GridspecPlot([s])
    p.render(
        save_plot=True,
        path_to_outfile=path_to_outfile
    )

def plot_3(ac, P_coag, P_frag):
    k = 0
    x = ac
    y1 = P_frag[k]
    y2 = P_coag[k]

    ax = plt.subplot()
    ax.set_xscale('log')
    # plt.scatter(x, y1)
    # plt.scatter(x, y2)
    # plt.plot(x, y1)
    # plt.plot(x, y2)
    plt.stackplot(
        x, [y2, y1], 
        colors=["red", "blue"],
        labels=["frag", "coag"],
    )

    plt.legend(loc="upper right")
    plt.xlim(ac[0], ac[-1])
    plt.ylim(0, 1)
    plt.show()
    plt.close()


for title, outcome_probability in outcome_probabilities.items():
    P_coag, P_frag = outcome_probability(cfg, dv)
    outcome_probabilities = {
        "P_coag": P_coag,
        "P_frag": P_frag,
    }

    plot_together(title, P_coag, P_frag)
    plot_separately(title, P_coag, P_frag)
    plot_3(ac, P_coag, P_frag)

