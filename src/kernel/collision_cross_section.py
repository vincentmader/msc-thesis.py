import numpy as np
from numpy import pi as PI

from disk.dust_particle import particle_radius_from_mass


def collision_cross_section(cfg, mg):
    if cfg.enable_physical_cross_sections is False:
        return np.ones(shape=[mg.N_x] * 2)

    masses = mg.grid_cell_centers()
    radii = particle_radius_from_mass(masses)

    sigma = np.ones(shape=[mg.N_x] * 2)
    for i, a_i in enumerate(radii):
        for j, a_j in enumerate(radii):
            sigma[i, j] = sigma_ij(a_i, a_j)
    return sigma


def sigma_ij(a_i, a_j):
    sigma_ij = PI * (a_i + a_j)**2
    return sigma_ij
