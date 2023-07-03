import numpy as np

from array_plotting import plot_array
from array_preparation import prepare_array
from array_flattening import flatten_array
from array_flattening import unflatten_array


ARR_N = 10
ARR_DIM = 3
ARR_SHAPE = [ARR_N]*ARR_DIM
ARR_MAX = ARR_N**2
CMAP_BOUNDS = 0, ARR_MAX


def prepare_p(arr):
    p = arr / np.sum(arr)
    return p


def main():

    arr = prepare_array(ARR_SHAPE)
    size = 10
    replace = False

    plot_array(arr, CMAP_BOUNDS)

    arr = flatten_array(arr)
    assert arr.all() == unflatten_array(flatten_array(arr), ARR_SHAPE).all()

    p = prepare_p(arr)

    a = np.random.choice(arr, size=size, replace=replace, p=p)
    print(a)

    arr = unflatten_array(arr, ARR_SHAPE)

    
if __name__ == "__main__":
    main()
