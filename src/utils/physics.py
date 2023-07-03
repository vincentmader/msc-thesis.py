import numpy as np

from constants import G, sigma_H2


def eddy_turnover_time(Omega_K):
    tau_ed = 1 / Omega_K
    return tau_ed


def reduced_mass(m_i, m_j):
    mu = (m_i * m_j) / (m_i + m_j)
    return mu


def mean_free_path(n):
    lambda_mfp = 1 / (n * sigma_H2)
    return lambda_mfp


def schmidt_nr(self, m):
    St = self.stokes_nr(m)
    Sc = 1 + St**2
    return Sc


def finite_difference(y, x):
    return (y[1:] - y[:-1]) / (x[1:] - x[:-1])


def kepler_frequency(r, M_star):
    Omega_K = np.sqrt(G * M_star / r**3)
    return Omega_K
