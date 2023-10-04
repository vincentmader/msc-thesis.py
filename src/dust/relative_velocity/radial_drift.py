import numpy as np

from utils import physics 


def dv_radial_drift(cfg, disk, disk_region):
    mc, ac  = disk.mg.bin_centers, disk.mg.particle_radii
    rho_s   = cfg.dust_particle_density
    rho_g   = disk_region.midplane_gas_volume_density
    u_th    = disk_region.thermal_velocity
    t_stop  = physics.stopping_time(rho_s, ac, rho_g, u_th)
    St      = physics.stokes_nr(ac, t_stop, rho_s)

    del_ln_P_g_del_ln_r      = disk_region.gas_pressure_gradient
    delr_Sigma_g_nu_g_sqrt_r = disk_region.delr_Sigma_g_nu_g_sqrt_r
    v_r = physics.u_r(
        cfg, disk_region, St, delr_Sigma_g_nu_g_sqrt_r, del_ln_P_g_del_ln_r
    )

    dv = np.zeros(shape=[disk.mg.N] * 2)
    for i, _ in enumerate(mc):
        v_i = v_r[i]
        for j, _ in enumerate(mc):
            v_j = v_r[j]

            dv[i, j] = np.abs(v_j - v_i)

    return dv


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
