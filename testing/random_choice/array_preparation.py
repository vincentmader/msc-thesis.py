import numpy as np


def prepare_array(shape):
    arr = np.zeros(shape=shape)
    for i in range(shape[0]):
        for j in range(shape[1]):
            for k in range(shape[2]):
                val = i * j * k
                arr[k, j, i] = val
    return arr
