from dust import particle_radius_from_mass
from utils.physics import mean_free_path, kepler_frequency


class DiskRegion:

    def __init__(self, cfg, disk):
        mg, rg = disk.mg, disk.rg
        mc, rc = mg.bin_centers, rg.bin_centers

        # Adjust distance from star, such that it's centered in its bin.
        r   = cfg.distance_to_star
        i_r = rg.index_from_value(r)
        r   = rc[i_r]

        # Load PPD properties as function of distance from star.
        Sigma_g    = disk.gas_surface_density
        M_star     = disk.stellar_mass
        T_mid      = disk.midplane_temperature
        c_s        = disk.sound_speed
        rho_g      = disk.midplane_gas_volume_density
        n          = disk.midplane_gas_volume_number_density
        u_th       = disk.thermal_velocity
        H_p        = disk.scale_height
        nu_mol     = disk.viscosity
        P_g_mid    = disk.midplane_gas_pressure
        rho_s      = cfg.dust_particle_density
        lambda_mfp = mean_free_path(n)
        Omega_K    = kepler_frequency(rc, M_star)
        v_K        = Omega_K * rc
        a          = particle_radius_from_mass(mc, rho_s)
        del_P_g_del_r            = disk.del_ln_P_g_del_ln_r
        delr_Sigma_g_nu_g_sqrt_r = disk.delr_Sigma_g_nu_g_sqrt_r

        # Determine PPD properties at radial distance `r`.
        self.mg, self.rg                 = mg, rg
        self.thermal_velocity            = u_th[i_r]
        self.midplane_gas_volume_density = rho_g[i_r]
        self.kepler_velocity             = v_K[i_r]
        self.gas_surface_density         = Sigma_g[i_r]
        self.viscosity                   = nu_mol[i_r]
        self.midplane_temperature        = T_mid[i_r]
        self.sound_speed                 = c_s[i_r]
        self.kepler_frequency            = Omega_K[i_r]
        self.kepler_velocity             = v_K[i_r]
        self.free_mean_path              = lambda_mfp[i_r]
        self.scale_height                = H_p[i_r]
        self.midplane_gas_pressure       = P_g_mid[i_r]
        self.distance_to_star            = r
        self.particle_radii              = a
        self.gas_pressure_gradient       = del_P_g_del_r[i_r]
        self.delr_Sigma_g_nu_g_sqrt_r    = delr_Sigma_g_nu_g_sqrt_r[i_r]
        # TODO ^ Rename?
