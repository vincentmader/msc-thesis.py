from typing import Optional

import numpy as np

from axis import DiscreteMassAxis, AxisLabelVariant, KernelErrorVariant
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
        R_coll: np.ndarray,
        axis_label_variant: Optional[AxisLabelVariant] = AxisLabelVariant.Radius,
        kernel_error_variant: Optional[KernelErrorVariant] = KernelErrorVariant.PercentPerCollision,
        *args, **kwargs
    ):
        err_matrix, err_total = test_mass_conservation(mg, K, R_coll, kernel_error_variant)

        if "title" not in kwargs.keys():
            kwargs["title"] = self.format_title(err_total, kernel_error_variant)

        super().__init__(cfg, mg, err_matrix, axis_label_variant=axis_label_variant, *args, **kwargs)

    def format_title(self, err_total, kernel_error_variant):
        if kernel_error_variant is KernelErrorVariant.KgPerSecond:
            dK_ij = r"$\sum_k m_k\cdot K_{kij}$"
            unit = "m$^3$ kg s$^{-1}$"
        elif kernel_error_variant is KernelErrorVariant.KgPerCollision:
            dK_ij = r"$\sum_k m_k\cdot\frac{K_{kij}}{R_{ij}}$"
            unit = "m$^3$ kg"
        elif kernel_error_variant is KernelErrorVariant.PercentPerCollision:
            dK_ij = r"$\sum_k \frac{m_k}{m_i+m_j}\cdot\frac{K_{kij}}{R_{ij}}$"
            unit = "m$^3$"
        else:
            raise Exception("NOTE: This can never happen.")

        dK = r"$\Delta K=\sqrt{\sum_{ij}\Delta K_{ij}^2}$ = " 
        title = r"$\Delta K_{ij}$ = " + f"{dK_ij}" + ", " + f"{dK} = {err_total:.2e} {unit}"
        return title

    def format_coord(self, x, y):
        try:
            if self.axis_variant is AxisLabelVariant.Radius:
                rho_s = self.cfg.dust_particle_density
                m_i = particle_mass_from_radius(y, rho_s)
                m_j = particle_mass_from_radius(x, rho_s)
                i = self.mg.index_from_value(m_i)
                j = self.mg.index_from_value(m_j)
            elif self.axis_variant is AxisLabelVariant.Radius:
                i = self.mg.index_from_value(y)
                j = self.mg.index_from_value(x)
            else:  # -> `AxisLabelVariant.Bin`
                i, j = int(y), int(x)
        except IndexError:
            return ""
        return f"{i = }, {j = }, sum_k K_kij = {self.z[i, j]:.2}"
