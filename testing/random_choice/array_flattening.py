import numpy as np


def flatten_array(arr):
    N = 1
    for i in arr.shape:
        N *= i

    shape = [N]
    return np.reshape(arr, shape)


def unflatten_array(arr, shape):
    return np.reshape(arr, shape)
