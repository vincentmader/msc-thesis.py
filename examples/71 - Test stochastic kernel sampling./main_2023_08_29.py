import os
import sys

import numpy as np
try:
    sys.path.append(os.path.join("..", "..", "src"))
    from config import Config
    from kernel import Kernel
except ModuleNotFoundError as e:
    raise e


def idx_2d_from_1d(index, N):
    i = index // N
    j = index % N
    return i, j

def sample_matrix_2d(
    matrix_2d: np.ndarray,
    prob_dist_2d: np.ndarray,
):
    # Assert 2D matrix.
    assert len(matrix_2d.shape) == 2, "matrix must be 2D"
    assert len(prob_dist_2d.shape) == 2, "matrix must be 2D"
    # Assert equal shapes for matrix & prob. dist.
    assert matrix_2d.shape[0] == prob_dist_2d.shape[0], "shape mismatch"
    assert matrix_2d.shape[1] == prob_dist_2d.shape[1], "shape mismatch"
    # Assert normalized prob. dist.
    assert prob_dist_2d.sum() - 1 < 1e-15, "prob_dist_2d must be normalized to 1"
    # Assert square matrices.
    N_i, N_j = matrix_2d.shape
    assert N_i == N_j

    # Create index matrix.
    indices_1d = np.arange(N**2)

    # Flatten: Convert from 2D to 1D.
    N_ij = N_i * N_j
    matrix_1d = matrix_2d.reshape([N_ij])
    prob_dist_1d = prob_dist_2d.reshape([N_ij])

    # Choose randomly from flattened matrix using flattened prob. dist.
    N_sample = N_i
    choices = np.random.choice(
        indices_1d,
        size=N_sample,
        replace=False,
        p=prob_dist_1d,
    )

    indices = [idx_2d_from_1d(idx, N_i) for idx in choices]
    return indices


def normalize(X):
    return X / np.sum(X)

def homogenous(N):
    P = np.ones(shape=(N, N))
    P = normalize(P)
    return P

cfg = Config()
kernel = Kernel(cfg)
K = kernel.K[25]

N = 50
i = np.arange(0, N, 1)
j = np.arange(0, N, 1)

# K = np.ones(shape=(N, N)) * i * j
P = homogenous(N)

indices = sample_matrix_2d(K, P)
print(indices)
