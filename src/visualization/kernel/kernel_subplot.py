from typing import Optional

import numpy as np

from axis import DiscreteMassAxis, KernelAxis
from visualization.base import PcolorMatrixSubplot


class KernelSubplot(PcolorMatrixSubplot):

    def __init__(
        self,
        mg: DiscreteMassAxis,
        K: np.ndarray,
        axis: Optional[KernelAxis] = KernelAxis.Radius,
        *args, **kwargs
    ):

        if axis is KernelAxis.Radius:
            kwargs["xlabel"] = kwargs["xlabel"]\
                if "xlabel" in kwargs.keys()\
                else "particle radius $a_j$ [m]"
            kwargs["ylabel"] = kwargs["ylabel"]\
                if "ylabel" in kwargs.keys()\
                else "particle radius $a_j$ [m]"
            ac = mg.particle_radii
            x, y = ac, ac
        elif axis is KernelAxis.Mass:
            kwargs["xlabel"] = kwargs["xlabel"]\
                if "xlabel" in kwargs.keys()\
                else "particle mass $a_j$ [kg]"
            kwargs["ylabel"] = kwargs["ylabel"]\
                if "ylabel" in kwargs.keys()\
                else "particle mass $a_i$ [kg]"
            mc = mg.bin_centers
            x, y = mc, mc
        else:  # -> `KernelAxis.Bin`
            kwargs["xlabel"] = kwargs["xlabel"]\
                if "xlabel" in kwargs.keys()\
                else "bin index $j$"
            kwargs["ylabel"] = kwargs["ylabel"]\
                if "ylabel" in kwargs.keys()\
                else "bin index $i$"
            i = np.linspace(0, mg.N, mg.N)
            x, y = i, i
        z = K

        super().__init__(x, y, z, *args, **kwargs)