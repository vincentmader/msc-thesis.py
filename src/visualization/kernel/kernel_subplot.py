from config import Config
from typing import Optional

import numpy as np

from axis import DiscreteMassAxis, AxisLabelVariant
from visualization.base import PcolorMatrixSubplot


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

        if axis_label_variant is AxisLabelVariant.Radius:
            kwargs["xlabel"] = kwargs["xlabel"]\
                if "xlabel" in kwargs.keys()\
                else "particle radius $a_j$ [m]"
            kwargs["ylabel"] = kwargs["ylabel"]\
                if "ylabel" in kwargs.keys()\
                else "particle radius $a_i$ [m]"
            ac = mg.particle_radii
            x, y = ac, ac
        elif axis_label_variant is AxisLabelVariant.Mass:
            kwargs["xlabel"] = kwargs["xlabel"]\
                if "xlabel" in kwargs.keys()\
                else "particle mass $a_j$ [kg]"
            kwargs["ylabel"] = kwargs["ylabel"]\
                if "ylabel" in kwargs.keys()\
                else "particle mass $a_i$ [kg]"
            mc = mg.bin_centers
            x, y = mc, mc
        else:  # -> `AxisLabelVariant.Bin`
            kwargs["xlabel"] = kwargs["xlabel"]\
                if "xlabel" in kwargs.keys()\
                else "bin index $j$"
            kwargs["ylabel"] = kwargs["ylabel"]\
                if "ylabel" in kwargs.keys()\
                else "bin index $i$"
            i = np.linspace(0, mg.N, mg.N)
            x, y = i, i

        self.cfg, self.mg, self.axis_variant = cfg, mg, axis_label_variant
        super().__init__(x, y, K, *args, **kwargs)
