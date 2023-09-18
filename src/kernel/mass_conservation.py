import numpy as np


def assert_cubic_kernel_shape(K):
    shape = K.shape
    assert len(shape) == 3
    assert shape[0] == shape[1]
    assert shape[1] == shape[2]
    return shape[0]


def test_mass_conservation(mg, K):
    K = np.array([0.5 * (K_k + K_k.T) for K_k in K])

    N_m = assert_cubic_kernel_shape(K)
    assert N_m == mg.N

    mc = mg.bin_centers
    dm = mg.bin_widths

    out = np.zeros(shape=[N_m] * 2)

    for i in range(N_m):
        for j in range(N_m):

            sum_ij = 0
            for k in range(N_m):
                m_k = mc[k]
                dm_k = dm[k]

                if mg.scale == "lin":
                    msg = f"dm_k != 1 on linear grid, is this on purpose? ({k=}, {dm_k=})"
                    assert np.abs(dm_k - 1) <= 1e-14, msg

                sum_ij += m_k * K[k, i, j]

            machine_precision_is_assured = abs(sum_ij)  # < 1e-12

            out[i, j] = machine_precision_is_assured

    return out
