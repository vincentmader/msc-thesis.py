from typing import Optional

import numpy as np

from axis import DiscreteMassAxis, KernelAxis
from config import Config
from dust import particle_mass_from_radius
from kernel.mass_conservation import test_mass_conservation
from visualization.kernel import KernelSubplot


class KernelMassConservationSubplot(KernelSubplot):
    __slots__ = ["K", "err_matrix", "err_total"]

    def __init__(
        self,
        cfg: Config,
        mg: DiscreteMassAxis,
        K: np.ndarray,
        axis_label_variant: Optional[KernelAxis] = KernelAxis.Radius, # TODO Rename? -> `KernelAxisLabelVariant`
        *args, **kwargs
    ):
        self.err_matrix, self.err_total = test_mass_conservation(mg, K)
        self.K = self.err_matrix

        kwargs["title"] = kwargs["title"]\
            if "title" in kwargs.keys()\
            else r"mass error $\Delta_{ij}=\sum_k m_k\cdot K_{kij}$, $\Delta=\sum_{ij}|\Delta_{ij}|$ = " + f"{self.err_total:.2e}" + " m$^3$s$^{-1}$kg"

        super().__init__(cfg, mg, self.K, axis_label_variant=axis_label_variant, *args, **kwargs)

    def format_coord(self, x, y):
        if self.axis_variant is KernelAxis.Radius:
            rho_s = self.cfg.dust_particle_density
            m_i = particle_mass_from_radius(y, rho_s)
            m_j = particle_mass_from_radius(x, rho_s)
            i = self.mg.index_from_value(m_i)
            j = self.mg.index_from_value(m_j)
        elif self.axis_variant is KernelAxis.Radius:
            i = self.mg.index_from_value(y)
            j = self.mg.index_from_value(x)
        else:  # -> `KernelAxis.Bin`
            i, j = int(y), int(x)

        return f"{i = }, {j = }, sum_k K_kij = {self.err_matrix[i, j]:.2}"
