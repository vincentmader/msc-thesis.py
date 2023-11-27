import numpy as np


def is_cubic(matrix: np.ndarray):
    shape = matrix.shape
    if len(shape) != 3:
        return False
    if shape[0] != shape[1]:
        return False
    if shape[1] != shape[2]:
        return False
    return True
