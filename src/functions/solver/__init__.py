import numpy as np


def fix_negative_densities(
    mc:     np.ndarray,
    N_dust: np.ndarray,
):
    idx_i = np.where(N_dust < 0)
    idx_k = np.where(N_dust >= 0)

    M_i = np.sum(N_dust[idx_i] * mc[idx_i])
    M_k = np.sum(N_dust[idx_k] * mc[idx_k])

    fac = np.abs(M_i) / np.abs(M_k)

    N_dust[idx_k] *= 1 - fac
    N_dust[idx_i] = 0

    return N_dust


def deriv(
    y: np.ndarray,
    x: np.ndarray,
):
    dx = x[1:] - x[:-1]
    dy = y[1:] - y[:-1]
    return dy / dx
