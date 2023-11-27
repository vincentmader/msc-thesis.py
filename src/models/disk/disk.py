from typing import Optional

from models.axis import DiscreteMassAxis, DiscreteRadialAxis
from config import Config
from functions.utils import physics


class Disk:

    def __init__(
        self,
        cfg:    Config,
        rg:     Optional[DiscreteRadialAxis]    = None,
        mg:     Optional[DiscreteMassAxis]      = None,
    ):

        self.mg = DiscreteMassAxis(cfg)     if mg is None else mg
        self.rg = DiscreteRadialAxis(cfg)   if rg is None else rg
        rc, rb = self.rg.bin_centers, self.rg.bin_boundaries

        L_star      = cfg.stellar_luminosity
        M_star      = cfg.stellar_mass
        M_disk      = M_star * cfg.disk_mass_ratio
        phi_fl      = cfg.flaring_angle
        T_mid       = physics.midplane_temperature(rc, L_star, phi_fl)
        c_s         = physics.sound_speed(rc, T_mid)
        H_p         = physics.scale_height(rc, M_star, c_s)
        Sigma_g     = physics.gas_surface_density(rb, M_disk)
        rho_g_mid   = physics.midplane_gas_volume_density(rc, Sigma_g, H_p)
        N_g_mid     = physics.midplane_gas_volume_number_density(rho_g_mid)
        P_g_mid     = physics.midplane_gas_pressure(rho_g_mid, c_s)
        u_th        = physics.thermal_velocity(c_s)
        lambda_mfp  = physics.mean_free_path(N_g_mid)
        nu          = physics.viscosity(u_th, lambda_mfp)
        del_ln_P_g_del_ln_r      = physics.del_ln_P_g_del_ln_r(rc, P_g_mid)
        delr_Sigma_g_nu_g_sqrt_r = physics.delr_Sigma_g_nu_g_sqrt_r(rc, Sigma_g, nu)

        self.dust_to_gas                        = cfg.dust_to_gas_ratio
        self.stellar_luminosity                 = L_star
        self.stellar_mass                       = M_star
        self.disk_mass                          = M_disk
        self.flaring_angle                      = phi_fl
        self.midplane_temperature               = T_mid
        self.sound_speed                        = c_s
        self.scale_height                       = H_p
        self.gas_surface_density                = Sigma_g
        self.midplane_gas_volume_density        = rho_g_mid
        self.midplane_gas_volume_number_density = N_g_mid
        self.midplane_gas_pressure              = P_g_mid
        self.thermal_velocity                   = u_th
        self.mean_free_path                     = lambda_mfp
        self.viscosity                          = nu
        self.del_ln_P_g_del_ln_r                = del_ln_P_g_del_ln_r
        self.delr_Sigma_g_nu_g_sqrt_r           = delr_Sigma_g_nu_g_sqrt_r

        self.mass_distribution = []  # TODO
