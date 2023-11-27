import numpy as np


def heaviside_theta(x):
    if x < 0:
        return 0
    if x > 0:
        return 1
    return 1 / 2


def root_mean_squared(arr):
    squares = [i**2 for i in arr]
    return sum(squares)**.5


def is_cubic(matrix: np.ndarray):
    shape = matrix.shape
    if len(shape) != 3:
        return False
    if shape[0] != shape[1]:
        return False
    if shape[1] != shape[2]:
        return False
    return True
