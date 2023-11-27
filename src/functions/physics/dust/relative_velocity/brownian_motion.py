import numpy as np
from numpy import pi as PI

from constants import k_B
from functions.utils.physics import reduced_mass


def dv_brownian_motion(cfg, disk, disk_region):
    mg = disk.mg
    mc = mg.bin_centers

    T_mid = disk_region.midplane_temperature
    dv = np.zeros(shape=[mg.N] * 2)
    for i, m_i in enumerate(mc):
        for j, m_j in enumerate(mc):

            beta = 1 / (k_B * T_mid)
            mu = reduced_mass(m_i, m_j)

            dv[i, j] = np.sqrt(8 / (PI * beta * mu))
            #        ^ adapted from 2010 Birnstiel & Dullemond

    return dv
