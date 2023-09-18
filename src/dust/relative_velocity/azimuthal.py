import numpy as np

from dust import particle_radius_from_mass

E_d = 0.5  # 2010 Birnstiel, TODO
E_d = 1    # Kees 2023-03-21


def dv_azimuthal(cfg, disk, disk_region):
    mg = disk.mg
    mc = mg.bin_centers
    rho_s = cfg.dust_particle_density
    ac = mg.particle_radii

    stopping_times = disk_region.stopping_time(ac, rho_s)
    stokes_nrs = disk_region.stokes_nr(ac, stopping_times, rho_s)
    u = u_n(disk_region)

    dv = np.zeros(shape=[mg.N] * 2)
    for i, _ in enumerate(mc):
        for j, _ in enumerate(mc):

            St_i = stokes_nrs[i]
            St_j = stokes_nrs[j]

            a = 1 / (St_i**2 + 1)
            b = 1 / (St_j**2 + 1)
            dv[i, j] = np.abs(u * (a - b))

    return dv


def u_n(disk_region):
    # Omega_K = disk_region.Omega_K
    # rho_g = disk_region.rho_g
    del_P_g_del_r = disk_region.gas_pressure_gradient

    # u_n = -E_d * del_P_g_del_r / (2 * rho_g * Omega_K)
    # ^ from 2010 Birnstiel, eq. 20

    c_s = disk_region.c_s
    v_K = disk_region.v_K
    u_n = - E_d / 2 * c_s**2 / v_K * del_P_g_del_r
    return u_n
