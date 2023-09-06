import os
import sys

import numpy as np
try:
    sys.path.append(os.path.join("..", "..", "src"))
    from config import Config
    from kernel import Kernel
except ModuleNotFoundError as e:
    raise e


def _weights(kernel: Kernel) -> np.ndarray:
    K = kernel.K
    mg = kernel.mg
    mc = mg.grid_cell_centers
    N_m = mg.N

    out = np.zeros(shape=(N_m, N_m))
    for k in range(N_m):
        out += mc[k] * np.abs(K[k])
    return out

def _probability(weights: np.ndarray) -> np.ndarray:
    S = weights.sum()
    weights = weights / S
    return weights

def _index_matrix(
    shape_2d: tuple[int, int],
) -> np.ndarray:
    N_i, N_j = shape_2d
    I = []
    for i in range(N_i):
        I_i = []
        for j in range(N_j):
            I_ij = i*N_i + j
            I_i.append(I_ij)
        I.append(I_i)
    return np.array(I)

def _sample_indices(
    P: np.ndarray,
    shape_2d: tuple[int, int],
) -> list[int]:
    I = _index_matrix(shape_2d)

    shape_1d = (shape_2d[0] * shape_2d[1])
    P = P.reshape(shape_1d)
    I = I.reshape(shape_1d)

    samples = np.random.choice(I, p=P, size=10)
    return list(samples)

def _ijs_from_indices(
    indices: list[int],
    shape_2d: tuple[int, int],
) -> list[tuple[int, int]]:
    
    out = []
    for idx in indices:
        i = idx // shape_2d[0]
        j = idx % shape_2d[0]
        out.append((i, j))
    return out

if __name__ == "__main__":
    cfg = Config()
    kernel = Kernel(cfg)

    W = _weights(kernel)
    P = _probability(W)

    # TODO Give `P` to integrator.
    # TODO In each time step, do `P -> P * N_i * N_j` or sth. like that.
    # TODO Call `_sample_indices` and `_ijs_from_indices`.
    # TODO Build kernel from only those `[(i,j), ...]`.
    # TODO Integrate with the new kernel.
    # TODO Compare results to earlier, as well as to analytical solutions.

    shape_2d = P.shape
    indices = _sample_indices(P, shape_2d)
    ijs = _ijs_from_indices(indices, shape_2d)
    print(ijs)
