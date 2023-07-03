import numpy as np

from array_plotting import plot_array
from array_preparation import prepare_array


ARR_N = 10
ARR_DIM = 3

ARR_SHAPE = [ARR_N]*ARR_DIM
FLAT_SHAPE = [ARR_N**ARR_DIM]

ARR_MAX = ARR_N**2
CMAP_BOUNDS = 0, ARR_MAX
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


def prepare_p(arr):
    p = arr / np.sum(arr)
    return p


def main():
    arr = prepare_array(ARR_SHAPE)
    plot_array(arr, CMAP_BOUNDS)

    arr = arr.reshape(FLAT_SHAPE)

    p = prepare_p(arr)
    choices = np.random.choice(arr, size=NR_OF_SAMPLES, replace=False, p=p)
    print(choices)

    index_tuples = index_tuple_matrix_from_3d_shape(ARR_SHAPE)
    shape = [ARR_SHAPE[0] * ARR_SHAPE[1] * ARR_SHAPE[2], 3]
    index_tuples = index_tuples.reshape(shape)

    choices = [index_tuples[int(c)] for c in choices]
    print(choices)


if __name__ == "__main__":
    main()
