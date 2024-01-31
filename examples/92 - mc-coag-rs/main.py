import numpy as np

import module


def _K_kij(ijs) -> np.ndarray:
    return "RUST: Build K_kij from pairs of (i, j)"

def _P_ij(N_k, W_ij) -> np.ndarray:
    return "RUST: Calculate P_ij from N_k and W_ij"
def _ijs(P_ij) -> list[tuple[int, int]]:
    return "RUST: Sample pairs of (i, j) using P_ij"

def _W_ij(R_ij) -> np.ndarray:
    K_kij  = _K_kij(None)
    KK_kij = R_ij * K_kij
    return "RUST: Calculate W_ij from KK_kij"

def _R_ij() -> np.ndarray:
    return "TODO: Calculate R_ij"
def _N_0() -> np.ndarray:
    return "TODO: Initialize N_k"

def step(N_k, R_ij, W_ij) -> np.ndarray:
    P_ij   = _P_ij(N_k, W_ij)
    ijs    = _ijs(P_ij)
    K_kij  = _K_kij(ijs)
    KK_kij = K_kij * R_ij
    return "TODO: Update N_k using KK_kij"

def run():
    R_ij   = _R_ij()
    W_ij   = _W_ij(R_ij)
    N_k    = _N_0()
    N_k    = step(N_k, R_ij, W_ij)
