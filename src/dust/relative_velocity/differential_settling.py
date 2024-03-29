import numpy as np


def dv_differential_settling(cfg, disk, disk_region):
    mg = disk.mg
    mc = mg.bin_centers

    rho_s = cfg.dust_particle_density

    Omega_K = disk_region.Omega_K
    stopping_times = disk_region.stopping_time(mc, rho_s)
    scale_heights = 1  # TODO
    settling_velocities = v_sett(stopping_times, Omega_K, scale_heights)

    dv = np.zeros(shape=[mg.N] * 2)
    for i, _ in enumerate(mc):
        for j, _ in enumerate(mc):

            v_i = settling_velocities[i]
            v_j = settling_velocities[j]
            dv[i, j] = np.abs(v_i - v_j)

    return np.zeros(shape=[mg.N] * 2)
    # ^ TODO Do not return this, but the actual `dv` instead.
    return dv


def v_sett(t_s, Omega_K, z):
    return t_s * Omega_K**2 * z

# def v1(m_i, m_j):  # From Kees' code in `dv_rel.py`
#     v_i = v_sett(m_i)
#     v_j = v_sett(m_j)
#     return np.abs(v_i - v_j)

# def v2(m_i, m_j):  # From 2007 Dullemond & Brauer
#     St_i = stokes_nr(m_i),
#     St_j = stokes_nr(m_j)
#     return z * Omega_K * np.abs(St_i / (1 + St_i) - St_j / (1 + St_j))

# def v3(): # From 2010 Birnstiel
#     h_i = 1  # todo
#     h_j = 1  # todo
#     v = Omega_K * np.abs(
#         h_i * min(St_i, 1 / 2) - h_j * min(St_j, 1 / 2)
#     )
