import numpy as np
from numpy import pi as PI

from constants import k_B
from utils.physics import reduced_mass


def dv_brownian_motion(cfg, disk, disk_region):
    mg = disk.mass_axis
    masses = mg.grid_cell_centers

    T_mid = disk_region.T_mid
    dv = np.zeros(shape=[mg.N] * 2)
    for i, m_i in enumerate(masses):
        for j, m_j in enumerate(masses):

            beta = 1 / (k_B * T_mid)
            mu = reduced_mass(m_i, m_j)

            dv[i, j] = np.sqrt(8 / (PI * beta * mu))
            #        ^ adapted from 2010 Birnstiel & Dullemond

    return dv
