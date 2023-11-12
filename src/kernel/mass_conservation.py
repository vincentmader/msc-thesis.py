import numpy as np

from utils.functions import is_cubic


def test_mass_conservation(mg, K, R_coll):  # TODO Turn this function into `Kernel` class method.

    # Assert correct shape of kernel matrix.
    assert is_cubic(K)
    N_m = K.shape[0]
    assert N_m == mg.N

    # Assert correct spacing on discretized mass axis.
    mc, dm = mg.bin_centers, mg.bin_widths
    if mg.scale == "lin":
        msg = f"dm_k != 1 on linear grid, is this on purpose?"
        assert np.abs(dm - 1).all() == 0, msg

    # Assert non-negative collision rates.
    assert R_coll.all() >= 0

    err_matrix = np.zeros(shape=[N_m] * 2)
    for i in range(N_m):
        for j in range(N_m):  # TODO Use `i+1` as max.?
            m_tot = mc[i] + mc[j]

            err_ij = 0
            for k in range(N_m):
                err_ij += mc[k] * K[k, i, j]
            err_matrix[i, j] = np.abs(err_ij / m_tot) / R_coll[i, j]
            #                TODO ^ Remove `abs()` here, plot +/- side by side.

    err_total = np.sum(err_matrix**2)**.5
    # err_total = np.sum(err_matrix)
    return err_matrix, err_total
