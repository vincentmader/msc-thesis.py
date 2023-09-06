import os
import sys

import numpy as np
try:
    sys.path.append(os.path.join("..", "..", "src"))
    from config import Config
    from disk import mass_distribution
    from kernel import Kernel
except ModuleNotFoundError as e:
    raise e


NR_OF_SAMPLES = 10


def index_tuple_matrix_from_3d_shape(shape):
    assert len(shape) == 3
    matrix = []
    for k in range(shape[0]):
        matrix_k = []
        for j in range(shape[1]):
            matrix_kj = []
            for i in range(shape[2]):
                matrix_kj.append((k, j, i))
            matrix_k.append(matrix_kj)
        matrix.append(matrix_k)
    return np.array(matrix)


def prepare_p(K, mg, n):
    mc = mg.grid_cell_centers
    K = np.abs(K)
    # Multiply kernel with particle abundances
    p = np.zeros(shape=K.shape)
    for k in range(mg.N):
        for i in range(mg.N):
            for j in range(mg.N):
                p[k, i, j] = K[k, i, j]  \
                    * n[i] * n[j] \
                    * mc[i] * mc[j]

    p = K / np.sum(K)
    return p


def main():
    # Define kernel.
    cfg = Config()
    kernel = Kernel(cfg)
    mg = kernel.mg
    K = kernel.K
    K = np.array([0.5 * (K_k + K_k.T) for K_k in K])

    # Initialize mass distribution.
    n0 = mass_distribution.dirac_delta(cfg)
    n = n0  # TODO

    # Create probability distribution from absolute kernel values.
    p = prepare_p(K, mg, n)

    # Flatten kernel for definition of probability distribution.
    FLAT_SHAPE = [K.shape[0] * K.shape[1] * K.shape[2]]
    # Flatten probability distribution
    p = p.reshape(FLAT_SHAPE)

    # Create 3d matrix of index tuples.
    index_tuples = index_tuple_matrix_from_3d_shape(K.shape)

    # Flatten index tuple matrix
    FLAT_SHAPE = [K.shape[0] * K.shape[1] * K.shape[2], 3]
    index_tuples = index_tuples.reshape(FLAT_SHAPE)

    # Choose randomly from tuple matrix.
    foo = range(len(index_tuples))  # <- TODO Rename.
    choices = np.random.choice(foo, size=NR_OF_SAMPLES, replace=False, p=p)
    choices = [index_tuples[int(c)] for c in choices]

    for c in choices:
        print(c)


if __name__ == "__main__":
    main()
