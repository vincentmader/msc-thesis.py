import numpy as np

from constants import DAYS_PER_YEAR


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


def format_time_as_years(
    time_in_seconds: float | int | np.float64,
) -> str:

    M = 60
    H = 60*M
    D = 24*H
    Y = DAYS_PER_YEAR * D
    KY = 1000 * Y
    MY = 1000 * KY

    if time_in_seconds >= MY:
        return f"{round(time_in_seconds / MY)} My"
    if time_in_seconds >= KY:
        return f"{round(time_in_seconds / KY)} ky"
    if time_in_seconds >= Y:
        return f"{round(time_in_seconds / Y)} y"
    if time_in_seconds >= D:
        return f"{round(time_in_seconds / D)} d"
    if time_in_seconds >= H:
        return f"{round(time_in_seconds / H)} h"
    if time_in_seconds >= M:
        return f"{round(time_in_seconds / M)} min"
    return f"{round(time_in_seconds)} s"
