from typing import Optional

import numpy as np

from axis import DiscreteMassAxis, KernelAxisLabelVariant
from config import Config
from dust import particle_mass_from_radius
from kernel.mass_conservation import test_mass_conservation
from visualization.kernel import KernelSubplot


class KernelMassConservationSubplot(KernelSubplot):
    __slots__ = ["err_matrix", "err_total"]

    def __init__(
        self,
        cfg: Config,
        mg: DiscreteMassAxis,
        K: np.ndarray,
        axis_label_variant: Optional[KernelAxisLabelVariant] = KernelAxisLabelVariant.Radius, # TODO Rename? -> `KernelAxisLabelVariant`
        *args, **kwargs
    ):
        err_matrix, err_total = test_mass_conservation(mg, K)

        title = r"error $\Delta K_{ij}=\sum_k \frac{m_k}{m_i+m_j}\cdot K_{kij}$, "
        title += r"$\Delta K=\sum_{ij}|\Delta K_{ij}|$ = " 
        title += f"{err_total:.2e}" + " m$^3$s$^{-1}$"

        kwargs["title"] = kwargs["title"] if "title" in kwargs.keys() else title
        super().__init__(cfg, mg, err_matrix, axis_label_variant=axis_label_variant, *args, **kwargs)

    def format_coord(self, x, y):
        if self.axis_variant is KernelAxisLabelVariant.Radius:
            rho_s = self.cfg.dust_particle_density
            m_i = particle_mass_from_radius(y, rho_s)
            m_j = particle_mass_from_radius(x, rho_s)
            i = self.mg.index_from_value(m_i)
            j = self.mg.index_from_value(m_j)
        elif self.axis_variant is KernelAxisLabelVariant.Radius:
            i = self.mg.index_from_value(y)
            j = self.mg.index_from_value(x)
        else:  # -> `KernelAxis.Bin`
            i, j = int(y), int(x)

        return f"{i = }, {j = }, sum_k K_kij = {self.z[i, j]:.2}"
