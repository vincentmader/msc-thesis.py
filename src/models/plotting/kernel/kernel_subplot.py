from config import Config
from typing import Optional

import numpy as np

from models.axis import DiscreteMassAxis, AxisLabelVariant
from functions.dust import particle_mass_from_radius
from models.plotting.base import PcolorMatrixSubplot


class KernelSubplot(PcolorMatrixSubplot):
    __slots__ = ["cfg", "mg", "axis_variant"]

    def __init__(
        self,
        cfg: Config,
        mg: DiscreteMassAxis,
        K: np.ndarray,
        axis_label_variant: Optional[AxisLabelVariant] = AxisLabelVariant.Radius,
        *args, **kwargs
    ):

        # TODO Fix x/y axes. Plot mass/bin/radius boundaries, not centers (?)
        if axis_label_variant is AxisLabelVariant.Radius:
            kwargs["xlabel"] = kwargs["xlabel"]\
                if "xlabel" in kwargs.keys()\
                else "Dust Particle Radius $a_j$ [m]"
            kwargs["ylabel"] = kwargs["ylabel"]\
                if "ylabel" in kwargs.keys()\
                else "Dust Particle Radius $a_i$ [m]"
            ac = mg.particle_radii
            x, y = ac, ac
        elif axis_label_variant is AxisLabelVariant.Mass:
            kwargs["xlabel"] = kwargs["xlabel"]\
                if "xlabel" in kwargs.keys()\
                else "Dust Particle Mass $m_j$ [kg]"
            kwargs["ylabel"] = kwargs["ylabel"]\
                if "ylabel" in kwargs.keys()\
                else "Dust Particle Mass $m_i$ [kg]"
            mc = mg.bin_centers
            x, y = mc, mc
        else:  # -> `AxisLabelVariant.Bin`
            kwargs["xlabel"] = kwargs["xlabel"]\
                if "xlabel" in kwargs.keys()\
                else "Bin Index $j$"
            kwargs["ylabel"] = kwargs["ylabel"]\
                if "ylabel" in kwargs.keys()\
                else "Bin Index $i$"
            i = np.linspace(0, mg.N, mg.N)
            x, y = i, i

        self.cfg, self.mg, self.axis_variant = cfg, mg, axis_label_variant
        super().__init__(x, y, K, *args, **kwargs)

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
        return f"k = {self.k}   {i = }   {j = }   K_kij = {self.z[self.k, i, j]:.2}"
