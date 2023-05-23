import numpy as np
from numpy import pi as PI

from config import Config
from constants import AU, sigma_SB, k_B, m_p
from utils.axis import DiscreteAxis
from utils.physics import kepler_frequency, mean_free_path, finite_difference
from .mass_grid import MassGrid


class Disk:
    def __init__(
        self,
        cfg: Config,
        radial_axis: DiscreteAxis,
        mass_axis: MassGrid,
    ):
        self.stellar_mass = cfg.stellar_mass
        self.stellar_luminosity = cfg.stellar_luminosity
        self.dust_to_gas = cfg.dust_to_gas_ratio
        self.flaring_angle = cfg.flaring_angle
        self.disk_mass = cfg.disk_mass_ratio * cfg.stellar_mass
        self.radial_axis = radial_axis
        self.mass_axis = mass_axis
        self.mass_distribution = []  # TODO

    def gas_surface_density(self, r):
        M_disk = self.disk_mass
        plsig = -1
        Sigma_g = (r[:-1] / AU)**plsig
        dsurf = PI * (r[1:]**2 - r[:-1]**2)
        mdummy = (Sigma_g * dsurf).sum()
        Sigma_g *= (M_disk / mdummy)
        # Sigma_g = np.append(Sigma_g, Sigma_g[-1])
        return Sigma_g

    def midplane_temperature(self, r):
        L_star = self.stellar_luminosity
        flang = self.flaring_angle  # (incidence angle of the light)
        T_mid = (flang * 0.5 * L_star / (4 * PI * r**2 * sigma_SB))**0.25
        return T_mid

    def sound_speed(self, r):
        T_mid = self.midplane_temperature(r=r)  # todo remove
        c_s = np.sqrt(k_B * T_mid / (2.3 * m_p))
        return c_s  # ^ from `notes_blackboard_talk.pdf`

    def scale_height(self, r):
        M_star = self.stellar_mass
        c_s = self.sound_speed(r=r)  # todo remove
        Omega_K = kepler_frequency(r, M_star)  # todo remove
        H_p = c_s / Omega_K
        return H_p

    def midplane_gas_volume_density(self, r, Sigma_g):
        H_p = self.scale_height(r)  # todo remove
        rho_g_mid = Sigma_g / (np.sqrt(2 * PI) * H_p)
        return rho_g_mid

    def midplane_gas_volume_number_density(self, rho_g):
        n = rho_g / (2.3 * m_p)
        return n

    def midplane_gas_pressure(self, r, Sigma_g):
        rho_g_mid = self.midplane_gas_volume_density(r, Sigma_g)  # todo remove
        c_s = self.sound_speed(r)  # todo remove
        P = rho_g_mid * c_s**2
        return P

    def thermal_velocity(self, r):
        c_s = self.sound_speed(r=r)  # todo remove
        return c_s * (8 / PI)**.5

    def viscosity(self, u_th, lambda_mfp):
        nu_mol = .5 * u_th * lambda_mfp
        return nu_mol

    def delr_Sigma_g_nu_g_sqrt_r(self, r, Sigma_g):
        u_th = self.thermal_velocity(r)
        rho_g = self.midplane_gas_volume_density(r, Sigma_g)
        n = self.midplane_gas_volume_number_density(rho_g)
        lambda_mfp = mean_free_path(n)
        nu_g = self.viscosity(u_th, lambda_mfp)

        Sigma_g_nu_g_sqrt_r = Sigma_g * nu_g * np.sqrt(r)
        delr_Sigma_g_nu_g_sqrt_r = finite_difference(Sigma_g_nu_g_sqrt_r, r)
        delr_Sigma_g_nu_g_sqrt_r = np.append(
            delr_Sigma_g_nu_g_sqrt_r, delr_Sigma_g_nu_g_sqrt_r[-1])  # TODO
        return delr_Sigma_g_nu_g_sqrt_r

    def del_ln_P_g_del_ln_r(self, r, Sigma_g):
        P_g = self.midplane_gas_pressure(r, Sigma_g)
        ln_P_g = np.log(P_g)
        ln_r = np.log(r)
        del_ln_P_g_del_ln_r = finite_difference(ln_P_g, ln_r)
        del_ln_P_g_del_ln_r = np.append(
            del_ln_P_g_del_ln_r, del_ln_P_g_del_ln_r[-1])  # TODO
        return del_ln_P_g_del_ln_r


def disk_mass_from_distribution(
    m: np.ndarray,
    n: np.ndarray,
):
    dm = m[1:] - m[:-1]
    m_tot = np.sum(m[:-1] * dm * n[:-1])  # TODO Is it correct like this?
    # m_tot = np.sum(m**2 * n)
    return m_tot


def disk_mass_error(
    m: np.ndarray,
    m_0: float
):
    return (m - m_0) / m_0
