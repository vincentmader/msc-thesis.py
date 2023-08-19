import os
import sys
try:
    sys.path.append(os.path.join("..", "..", "src"))
    from axis import DiscreteMassAxis, DiscreteRadialAxis
    from collision import collision_outcome_probabilities_from_maxwell_boltzmann 
    from collision import collision_outcome_probabilities_from_cutoff_velocity
    from config import Config
    from disk import Disk, DiskRegion
    from dust import particle_radius_from_mass
    from dust.relative_velocity import relative_velocity
    from visualization.kernel_v3_2023_08_14.pcolor_matrix_subplot import PcolorMatrixSubplot
    from visualization.kernel_v3_2023_08_14.gridspec_plot import GridspecPlot
except ModuleNotFoundError as e:
    raise e


# Define kernel configuration.
cfg = Config()
rho_s = cfg.dust_particle_density

# Define discrete axis for radial distance from star, as well as for mass.
rg = DiscreteRadialAxis(cfg)
mg = DiscreteMassAxis(cfg)
mc = mg.grid_cell_centers
ac = particle_radius_from_mass(mc, rho_s)

# Define disk, the position of interest in it, & the disk properties there.
disk = Disk(cfg, rg, mg)
disk_region = DiskRegion(cfg, disk)

dv = relative_velocity(cfg, disk, disk_region)

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

    s1 = PcolorMatrixSubplot(
        ac, ac, P_coag,
        title="dust particle coagulation probability $P_{coag}$",
        xlabel="particle radius $a_j$",
        ylabel="particle radius $a_i$",
        scales=("log", "log", "lin"),
    )
    s2 = PcolorMatrixSubplot(
        ac, ac, P_frag,
        title="dust particle fragmentation probability $P_{frag}$",
        xlabel="particle radius $a_j$",
        ylabel="particle radius $a_i$",
        scales=("log", "log", "lin"),
    )
    
    subplots = [s1, s2]
    
    p = GridspecPlot(subplots)
    p.render()
