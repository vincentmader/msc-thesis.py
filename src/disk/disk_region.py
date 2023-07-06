import numpy as np
from numpy import pi as PI

from constants import rho_s
from dust import particle_radius_from_mass
from utils.physics import mean_free_path, eddy_turnover_time, finite_difference, kepler_frequency


class DiskRegion:
    def __init__(self, cfg, disk):

        mg = disk.mass_axis
        masses = mg.grid_cell_centers  # TODO

        rg = disk.radial_axis
        rc = rg.grid_cell_centers
        r = cfg.distance_to_star
        i_r = rg.index_from_value(r)  # TODO
        r = rc[i_r]

        Sigma_g = disk.gas_surface_density
        radii = particle_radius_from_mass(masses)

        M_star = disk.stellar_mass
        T_mid = disk.midplane_temperature
        c_s = disk.sound_speed
        rho_g = disk.midplane_gas_volume_density
        Omega_K = kepler_frequency(rc, M_star)
        v_K = Omega_K * rc
        n = disk.midplane_gas_volume_number_density
        lambda_mfp = mean_free_path(n)
        u_th = disk.thermal_velocity
        H_p = disk.scale_height
        nu_mol = disk.viscosity

        P = disk.midplane_gas_pressure
        logP, logr = np.log(P), np.log(rc)
        gas_pressure_gradient = finite_difference(logP, logr)

        delr_Sigma_g_nu_g_sqrt_r = disk.delr_Sigma_g_nu_g_sqrt_r[i_r]

        self.r = r
        self.T_mid = T_mid[i_r]
        self.v_K = v_K[i_r]
        self.c_s = c_s[i_r]
        self.Sigma_g = Sigma_g[i_r]  # TODO index correct?
        self.Omega_K = Omega_K[i_r]
        self.lambda_mfp = lambda_mfp[i_r]
        self.u_th = u_th[i_r]
        self.H_p = H_p[i_r]
        self.P = P[i_r]
        self.nu_mol = nu_mol[i_r]
        self.rho_g = rho_g[i_r]
        self.radii = radii

        self.gas_pressure_gradient = gas_pressure_gradient[i_r]
        # ^ TODO Rename (log/log deriv)
        self.delr_Sigma_g_nu_g_sqrt_r = delr_Sigma_g_nu_g_sqrt_r
        # TODO Rename?

    def stopping_time(self, radii):
        u_th = self.u_th
        rho_g = self.rho_g

        t_s = (rho_s * radii) / (rho_g * u_th)

        # TODO Implement Stokes regimes.
        # Re = self.reynolds_nr(radii, u)
        # nu_mol = self.nu_mol
        # res = np.zeros(shape=(len(radii)))
        # u = 1  # TODO
        #     if Re < 1:
        #         t_s_i = (2 * rho_s * a**2) / (9 * nu_mol * rho_g)
        #     if 1 < Re < 800:
        #         t_s_i = (2**0.6 * rho_s * a**1.6) /
        #         (9 * nu_mol**0.6 * rho_g**1.4 * u**0.4)
        #     if Re > 800:
        #         return (6 * rho_s * a) / (rho_g * u)
        #     # todo: when epstein? see 2010 Birnstiel
        return t_s

    def stokes_nr(self, a, t_stop):
        Sigma_g = self.Sigma_g
        Omega_K = self.Omega_K
        tau_ed = eddy_turnover_time(Omega_K)
        lambda_mfp = self.lambda_mfp

        St_v1 = rho_s * a / Sigma_g * PI / 2    # 2010 Birnstiel
        St_v2 = t_stop / tau_ed                 # 2010 Birnstiel

        # TODO Distinguish cases.
        if True:  # a < 9 / 4 * lambda_mfp:     # 2010 Birnstiel
            return St_v1
        else:
            return St_v2

    def reynolds_nr(self, a, u):  # TODO Don't use Kepler here!
        nu_mol = self.nu_mol
        Re = 2 * a * u / nu_mol
        return Re
