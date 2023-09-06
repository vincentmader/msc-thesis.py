from typing import Optional

import numpy as np
from numpy import pi as PI

from axis import DiscreteMassAxis, DiscreteRadialAxis
from config import Config
from constants import AU, sigma_SB, k_B, m_p
from utils.physics import kepler_frequency, mean_free_path, finite_difference


class Disk:

    def __init__(
        self,
        cfg: Config,
        radial_axis: Optional[DiscreteRadialAxis] = None,
        mass_axis: Optional[DiscreteMassAxis] = None,
    ):
        self.mass_axis = DiscreteMassAxis(cfg) if mass_axis is None else mass_axis
        self.radial_axis = DiscreteRadialAxis(cfg) if radial_axis is None else radial_axis

        self.stellar_luminosity = cfg.stellar_luminosity
        self.stellar_mass = cfg.stellar_mass
        self.disk_mass = cfg.disk_mass_ratio * cfg.stellar_mass
        self.dust_to_gas = cfg.dust_to_gas_ratio

        self.flaring_angle = cfg.flaring_angle
        self.mass_distribution = []  # TODO

        self.gas_surface_density= self._gas_surface_density()
        self.midplane_temperature = self._midplane_temperature()
        self.sound_speed = self._sound_speed()
        self.scale_height = self._scale_height()
        self.midplane_gas_volume_density = self._midplane_gas_volume_density()
        self.midplane_gas_volume_number_density = self._midplane_gas_volume_number_density()
        self.midplane_gas_pressure = self._midplane_gas_pressure()
        self.thermal_velocity = self._thermal_velocity()
        self.viscosity = self._viscosity()
        self.delr_Sigma_g_nu_g_sqrt_r = self._delr_Sigma_g_nu_g_sqrt_r()
        self.del_ln_P_g_del_ln_r = self._del_ln_P_g_del_ln_r()

    def _gas_surface_density(
        self, 
        r: Optional[np.ndarray] = None,
        M_disk: Optional[np.float64] = None,
    ):
        r = self.radial_axis.grid_cell_boundaries if r is None else r # NOTE: bounds!
        M_disk = self.disk_mass if M_disk is None else M_disk

        plsig = -1
        Sigma_g = (r[:-1] / AU)**plsig
        dsurf = PI * (r[1:]**2 - r[:-1]**2)
        mdummy = (Sigma_g * dsurf).sum()
        Sigma_g *= (M_disk / mdummy)
        # Sigma_g = np.append(Sigma_g, Sigma_g[-1])
        return Sigma_g

    def _midplane_temperature(
        self, 
        r: Optional[np.ndarray] = None,
        L_star: Optional[np.float64] = None,
        flang: Optional[np.float64] = None,
    ):
        r = self.radial_axis.grid_cell_centers if r is None else r
        L_star = self.stellar_luminosity if L_star is None else L_star
        flang = self.flaring_angle if flang is None else flang
        # ^ (incidence angle of the light)

        T_mid = (flang * 0.5 * L_star / (4 * PI * r**2 * sigma_SB))**0.25
        return T_mid

    def _sound_speed(
        self, 
        r: Optional[np.ndarray] = None,
        T_mid: Optional[np.ndarray] = None,
    ):
        r = self.radial_axis.grid_cell_centers if r is None else r
        T_mid = self.midplane_temperature if T_mid is None else T_mid

        c_s = np.sqrt(k_B * T_mid / (2.3 * m_p))
        return c_s  # ^ from `notes_blackboard_talk.pdf`

    def _scale_height(
        self, 
        r: Optional[np.ndarray] = None,
        M_star: Optional[np.float64] = None,
        c_s: Optional[np.ndarray] = None
    ):
        r = self.radial_axis.grid_cell_centers if r is None else r
        M_star = self.stellar_mass if M_star is None else M_star
        c_s = self.sound_speed if c_s is None else c_s

        Omega_K = kepler_frequency(r, M_star)  # todo remove
        H_p = c_s / Omega_K
        return H_p

    def _midplane_gas_volume_density(
        self, 
        r: Optional[np.ndarray] = None,
        Sigma_g: Optional[np.ndarray] = None,
        H_p: Optional[np.ndarray] = None
    ):
        r = self.radial_axis.grid_cell_centers if r is None else r
        Sigma_g = self.gas_surface_density if Sigma_g is None else Sigma_g
        H_p = self.scale_height if H_p is None else H_p

        rho_g_mid = Sigma_g / (np.sqrt(2 * PI) * H_p)
        return rho_g_mid

    def _midplane_gas_volume_number_density(
        self, 
        rho_g: Optional[np.ndarray] = None,
    ):
        rho_g = self.midplane_gas_volume_density if rho_g is None else rho_g

        n = rho_g / (2.3 * m_p)
        return n

    def _midplane_gas_pressure(
        self, 
        rho_g_mid: Optional[np.ndarray] = None,
        c_s: Optional[np.ndarray] = None,
    ):
        rho_g_mid = self.midplane_gas_volume_density if rho_g_mid is None else rho_g_mid
        c_s = self.sound_speed if c_s is None else c_s

        P = rho_g_mid * c_s**2
        return P

    def _thermal_velocity(
        self, 
        c_s: Optional[np.ndarray] = None,
    ):
        c_s = self.sound_speed if c_s is None else c_s

        return c_s * (8 / PI)**.5

    def _viscosity(
        self, 
        u_th: Optional[np.ndarray] = None,  
        lambda_mfp: Optional[np.ndarray] = None, 
    ):
        u_th = self.thermal_velocity if u_th is None else u_th
        if lambda_mfp is None:
            n = self.midplane_gas_volume_number_density
            lambda_mfp = mean_free_path(n)  # TODO mean free path of what?

        nu_mol = .5 * u_th * lambda_mfp
        return nu_mol

    def _delr_Sigma_g_nu_g_sqrt_r(
        self, 
        r: Optional[np.ndarray] = None,
        Sigma_g: Optional[np.ndarray] = None,
        nu_g: Optional[np.ndarray] = None,
    ):
        r = self.radial_axis.grid_cell_centers if r is None else r
        Sigma_g = self.gas_surface_density if Sigma_g is None else Sigma_g
        nu_g = self.viscosity if nu_g is None else nu_g

        Sigma_g_nu_g_sqrt_r = Sigma_g * nu_g * np.sqrt(r)

        delr_Sigma_g_nu_g_sqrt_r = finite_difference(Sigma_g_nu_g_sqrt_r, r)
        delr_Sigma_g_nu_g_sqrt_r = np.append(
            delr_Sigma_g_nu_g_sqrt_r, delr_Sigma_g_nu_g_sqrt_r[-1]
        )  # TODO
        return delr_Sigma_g_nu_g_sqrt_r

    def _del_ln_P_g_del_ln_r(
        self, 
        r: Optional[np.ndarray] = None,
        P_g: Optional[np.ndarray] = None,
    ):
        r = self.radial_axis.grid_cell_centers if r is None else r
        P_g = self.midplane_gas_pressure if P_g is None else P_g

        ln_P_g = np.log(P_g)
        ln_r = np.log(r)

        del_ln_P_g_del_ln_r = finite_difference(ln_P_g, ln_r)

        del_ln_P_g_del_ln_r = np.append(
            del_ln_P_g_del_ln_r, del_ln_P_g_del_ln_r[-1]
        )  # TODO
        return del_ln_P_g_del_ln_r


def disk_mass_from_distribution(
    n: np.ndarray,
    m: np.ndarray,
    dm: np.ndarray,
):
    m_tot = np.sum(m * dm * n)  # TODO Is it correct like this?
    return m_tot


def disk_mass_error(
    m: np.ndarray,
    m_0: float
):
    return (m - m_0) / m_0
