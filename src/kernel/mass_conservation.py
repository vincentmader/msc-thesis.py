import numpy as np


def assert_cubic_kernel_shape(K):
    shape = K.shape
    assert len(shape) == 3
    assert shape[0] == shape[1]
    assert shape[1] == shape[2]
    return shape[0]


def test_mass_conservation(mg, K):  # TODO Turn this function into `Kernel` class method.
    K = np.array([0.5 * (K_k + K_k.T) for K_k in K])  
    # ^ TODO Why such different results when this line is commented out?

    N_m = assert_cubic_kernel_shape(K)
    assert N_m == mg.N

    mc, dm = mg.bin_centers, mg.bin_widths
    if mg.scale == "lin":
        msg = f"dm_k != 1 on linear grid, is this on purpose?"
        assert np.abs(dm - 1).all() == 0, msg

    err_matrix = np.zeros(shape=[N_m] * 2)
    for i in range(N_m):
        for j in range(N_m):  # TODO Use `i+1` as max.?
            m_tot = mc[i] + mc[j]

            err_ij = 0
            for k in range(N_m):
                err_ij += mc[k] * K[k, i, j]
            err_matrix[i, j] = np.abs(err_ij / m_tot)  
            #                TODO ^ Remove `abs()` here, plot +/- side by side.

    err_total = np.sum(err_matrix**2)**.5
    # err_total = np.sum(err_matrix)
    return err_matrix, err_total
