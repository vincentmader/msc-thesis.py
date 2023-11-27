import numpy as np

from constants import DAYS_PER_YEAR


def format_seconds_as_years(
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
