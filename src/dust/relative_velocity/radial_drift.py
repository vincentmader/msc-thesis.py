import numpy as np

from disk.dust_particle import particle_radius_from_mass


def dv_radial_drift(cfg, disk, disk_region):
    mg = disk.mass_axis

    del_ln_P_g_del_ln_r = disk_region.gas_pressure_gradient
    delr_Sigma_g_nu_g_sqrt_r = disk_region.delr_Sigma_g_nu_g_sqrt_r

    mc = mg.grid_cell_centers  # TODO Use bounds or centers?
    radii = particle_radius_from_mass(mc)
    stopping_times = disk_region.stopping_time(radii)
    stokes_nrs = disk_region.stokes_nr(radii, stopping_times)

    v_r = u_r(cfg, disk_region, stokes_nrs,
              delr_Sigma_g_nu_g_sqrt_r, del_ln_P_g_del_ln_r)

    dv = np.zeros(shape=[mg.N] * 2)
    for i, _ in enumerate(mc):
        v_i = v_r[i]
        for j, _ in enumerate(mc):
            v_j = v_r[j]

            dv[i, j] = np.abs(v_j - v_i)

    return dv


def u_r(cfg, disk_region, stokes_nrs, delr_Sigma_g_nu_g_sqrt_r, del_ln_P_g_del_ln_r):
    """Radial Dust Velocity"""
    St = stokes_nrs
    c_s = disk_region.c_s
    v_K = disk_region.v_K
    Sigma_g = disk_region.Sigma_g
    r = cfg.distance_to_star

    v_g = u_g(r, Sigma_g, delr_Sigma_g_nu_g_sqrt_r)
    v_n = u_n(c_s, v_K, del_ln_P_g_del_ln_r)
    v_r = (v_g - 2 * St * v_n) / (1 + St**2)
    return v_r


def u_g(r, Sigma_g, delr_Sigma_g_nu_g_sqrt_r):
    """Radial Gas Velocity"""
    u_g = -3 / (Sigma_g * np.sqrt(r)) * delr_Sigma_g_nu_g_sqrt_r
    return u_g


def u_n(c_s, v_K, del_ln_P_g_del_ln_r):
    """Maximum drift velocity of a particle"""
    E_d = 1
    u_n = - E_d / 2 * c_s**2 / v_K * del_ln_P_g_del_ln_r
    # u_n = -del_P_g_del_r * E_d / (2 * rho_g * Omega_K)
    return u_n


# def v_r(disk_region, masses, del_P_g_del_r):
#     radii = particle_radius_from_mass(masses)
#     stopping_times = disk_region.stopping_time(radii)

#     # TODO Handle both small & large particles.
#     v = v_r_large(disk_region, del_P_g_del_r, stopping_times)
#     return v


# def v_r_small(disk_region, del_P_g_del_r, stopping_times):
#     c_s = disk_region.c_s
#     r = disk_region.r

#     v = del_P_g_del_r * c_s**2 / r * stopping_times
#     return v


# def v_r_large(disk_region, del_P_g_del_r, stopping_times):
#     c_s = disk_region.c_s
#     r = disk_region.r
#     Omega_K = disk_region.Omega_K
#     v_K = Omega_K * r

#     v = del_P_g_del_r * c_s**2 / v_K**2 * r / stopping_times
#     return v


# # def u_r(m):
# #     St = stokes_nr(m)
# #     return u_g(R) / (1 + St**2) - 2 * u_n() / (St + 1 / St)

# # def u_g(r):
# #     return -3 / (Sigma_g * np.sqrt(r)) * delr_Sigma_g_nu_g_sqrt_r

# # def u_n():
# #     return -E_d * dPdr / (2 * rho_g * Omega_K)
