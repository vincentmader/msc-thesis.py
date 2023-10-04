import numpy as np
from numpy import pi as PI

from constants import G, AU, sigma_SB, sigma_H2, k_B, m_p


def eddy_turnover_time(Omega_K):
    tau_ed = 1 / Omega_K
    return tau_ed


def reduced_mass(m_i, m_j):
    mu = (m_i * m_j) / (m_i + m_j)
    return mu


def mean_free_path(n):
    lambda_mfp = 1 / (n * sigma_H2)
    return lambda_mfp


# def schmidt_nr(self, m):
#     St = self.stokes_nr(m)
#     Sc = 1 + St**2
#     return Sc


def finite_difference(y, x):
    return (y[1:] - y[:-1]) / (x[1:] - x[:-1])


def kepler_frequency(r, M_star):
    Omega_K = np.sqrt(G * M_star / r**3)
    return Omega_K


def disk_mass_from_distribution(n: np.ndarray, m: np.ndarray, dm: np.ndarray):
    m_tot = np.sum(m * dm * n)  # TODO Is it correct like this?
    return m_tot


def disk_mass_error(m: np.ndarray, m_0: float):
    return (m - m_0) / m_0


def thermal_velocity(c_s: np.ndarray):
    u_th = c_s * (8 / PI)**.5
    return u_th


def midplane_gas_pressure(rho_g_mid: np.ndarray, c_s: np.ndarray):
    P = rho_g_mid * c_s**2
    return P


def midplane_gas_volume_number_density(rho_g_mid: np.ndarray):
    n = rho_g_mid / (2.3 * m_p)
    return n


def sound_speed(r: np.ndarray, T_mid: np.ndarray):
    c_s = np.sqrt(k_B * T_mid / (2.3 * m_p))
    return c_s  # ^ from `notes_blackboard_talk.pdf`


def scale_height(r: np.ndarray, M_star: float, c_s: np.ndarray):
    Omega_K = kepler_frequency(r, M_star)  # TODO Remove?
    H_p = c_s / Omega_K
    return H_p


def gas_surface_density(r: np.ndarray, M_disk: float):
    plsig = -1
    Sigma_g = (r[:-1] / AU)**plsig
    dsurf = PI * (r[1:]**2 - r[:-1]**2)
    mdummy = (Sigma_g * dsurf).sum()
    Sigma_g *= (M_disk / mdummy)
    # Sigma_g = np.append(Sigma_g, Sigma_g[-1])
    return Sigma_g


def midplane_temperature(r: np.ndarray, L_star: float, flang: float):
    T_mid = (flang * 0.5 * L_star / (4 * PI * r**2 * sigma_SB))**0.25
    return T_mid  # ^ (flang = incidence angle of the light)


def viscosity(u_th: np.ndarray, lambda_mfp: np.ndarray):
    nu_mol = .5 * u_th * lambda_mfp
    return nu_mol


def midplane_gas_volume_density(r: np.ndarray, Sigma_g: np.ndarray, H_p: np.ndarray):
    rho_g_mid = Sigma_g / (np.sqrt(2 * PI) * H_p)
    return rho_g_mid


def delr_Sigma_g_nu_g_sqrt_r(r: np.ndarray, Sigma_g: np.ndarray, nu_g: np.ndarray):
    Sigma_g_nu_g_sqrt_r = Sigma_g * nu_g * np.sqrt(r)
    delr_Sigma_g_nu_g_sqrt_r = finite_difference(Sigma_g_nu_g_sqrt_r, r)
    # TODO Make sure the output array has correct shape (append one more?).
    return delr_Sigma_g_nu_g_sqrt_r


def del_ln_P_g_del_ln_r(r: np.ndarray, P_g: np.ndarray):
    ln_P_g = np.log(P_g)
    ln_r = np.log(r)
    del_ln_P_g_del_ln_r = finite_difference(ln_P_g, ln_r)
    # TODO Make sure the output array has correct shape (append one more?).
    return del_ln_P_g_del_ln_r


def stopping_time(rho_s, a, rho_g, u_th):
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
    t_s = (rho_s * a) / (rho_g * u_th)
    return t_s


def stokes_nr(rho_s, a, Sigma_g):
    # TODO Distinguish cases.
    if True:  # a < 9 / 4 * lambda_mfp:      # 2010 Birnstiel
        St = rho_s * a / Sigma_g * PI / 2    # 2010 Birnstiel
        return St
    else:
        St = t_stop / tau_ed                 # 2010 Birnstiel
        return St


def reynolds_nr(a, u, nu):  # NOTE Don't use Kepler here!
    Re = 2 * a * u / nu
    return Re
