import numpy as np

from utils import physics


E_d = 0.5  # 2010 Birnstiel, TODO
E_d = 1    # Kees 2023-03-21


def dv_azimuthal(cfg, disk, disk_region):
    mc, ac  = disk.mg.bin_centers, disk.mg.particle_radii
    rho_s   = cfg.dust_particle_density
    rho_g   = disk_region.midplane_gas_volume_density
    u_th    = disk_region.thermal_velocity
    t_stop  = physics.stopping_time(rho_s, ac, rho_g, u_th)
    St      = physics.stokes_nr(ac, t_stop, rho_s)

    c_s     = disk_region.sound_speed
    v_K     = disk_region.kepler_velocity
    del_ln_P_g_del_ln_r = disk_region.gas_pressure_gradient
    u_n     = physics.u_n(c_s, v_K, del_ln_P_g_del_ln_r)

    dv = np.zeros(shape=[disk.mg.N] * 2)
    for i, _ in enumerate(mc):
        for j, _ in enumerate(mc):

            St_i, St_j = St[i], St[j]
            a = 1 / (St_i**2 + 1)
            b = 1 / (St_j**2 + 1)
            dv[i, j] = np.abs(u_n * (a - b))

    return dv
