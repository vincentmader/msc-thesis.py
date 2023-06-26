import numpy as np


def assert_cubic_kernel_shape(K):
    shape = K.shape
    assert len(shape) == 3
    assert shape[0] == shape[1]
    assert shape[1] == shape[2]
    return shape[0]


def test_mass_conservation(kernel):
    K = kernel.K
    mg = kernel.mg

    K = np.array([0.5 * (K_k + K_k.T) for K_k in K])

    N_m = assert_cubic_kernel_shape(K)
    assert N_m == mg.N

    mc = mg.grid_cell_centers
    dm = mg.grid_cell_widths

    out = np.zeros(shape=[N_m] * 2)

    for i in range(N_m):
        for j in range(N_m):

            sum_ij = 0
            for k in range(N_m):
                m_k = mc[k]
                dm_k = dm[k]

                if kernel.cfg.mass_axis_scale == "lin" and np.abs(dm_k - 1) > 1e-14:
                    msg = f"Grid spacing dm_k != 1 on linear grid, is this on purpose? (for {k=} -> {dm_k=})"
                    raise Exception(msg)

                sum_ij += m_k * dm_k * K[k, i, j]

            machine_precision_is_assured = abs(sum_ij)  # < 1e-12

            out[i, j] = machine_precision_is_assured

    return out
