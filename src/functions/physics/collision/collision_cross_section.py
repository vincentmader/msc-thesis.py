import numpy as np
from numpy import pi as PI

from functions.physics.dust import particle_radius_from_mass


def collision_cross_section(cfg, mg):
    if cfg.enable_physical_collisions is False:
        return np.ones(shape=[mg.N] * 2)

    rho_s = cfg.dust_particle_density
    mc = mg.bin_centers
    radii = particle_radius_from_mass(mc, rho_s)

    sigma = np.ones(shape=[mg.N] * 2)
    for i, a_i in enumerate(radii):
        for j, a_j in enumerate(radii):
            sigma[i, j] = sigma_ij(a_i, a_j)
    return sigma


def sigma_ij(a_i, a_j):
    sigma_ij = PI * (a_i + a_j)**2
    return sigma_ij
