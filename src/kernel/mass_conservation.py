from typing import Optional

import numpy as np

from axis import KernelErrorVariant
from utils.functions import is_cubic


def test_mass_conservation(
    mg, K, R_coll, 
    kernel_error_variant: Optional[KernelErrorVariant] = KernelErrorVariant.PercentPerCollision,
):  
    # TODO Turn this function into `Kernel` class method.

    # PERFORM ASSERTIONS.
    # ─────────────────────────────────────────────────────────────────────────

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

    # CALCULATE KERNEL MASS ERROR FOR ALL COLLISIONS `(i, j)`.
    # ─────────────────────────────────────────────────────────────────────────

    # Initialize 2d matrix for error in collision `(i, j)`.
    E = np.zeros(shape=[N_m] * 2)

    # Calculate error for all collisions `(i, j)`.
    for i in range(N_m):
        for j in range(N_m):  # <- TODO Use `i+1` as max.?
            m_tot = mc[i] + mc[j]

            E_ij = 0
            for k in range(N_m):
                E_ij += mc[k] * K[k, i, j]
            E[i, j] = np.abs(E_ij)
            #  TODO ^ Remove `abs()` here, plot +/- side by side.
            if kernel_error_variant == KernelErrorVariant.KgPerCollision:
                E[i, j] /= R_coll[i, j]
            if kernel_error_variant == KernelErrorVariant.PercentPerCollision:
                E[i, j] /= R_coll[i, j] * m_tot / 100

    # Calculate total error.
    E_tot = np.sum(E**2)**.5
    return E, E_tot
