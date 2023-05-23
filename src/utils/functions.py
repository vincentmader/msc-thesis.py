def heaviside_theta(x):
    if x < 0:
        return 0
    if x > 0:
        return 1
    return 1 / 2


def root_mean_squared(arr):
    # TODO Use np.sum instead?
    return sum([i**2 for i in arr])**.5
