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
            xlabel = "particle radius $a_j$"
            ylabel = "particle radius $a_i$"
            ac = mg.particle_radii
            x, y = ac, ac
        elif axis is KernelAxis.Mass:
            xlabel = "particle mass $a_j$"
            ylabel = "particle mass $a_i$"
            mc = mg.grid_cell_centers
            x, y = mc, mc
        else:  # -> `KernelAxis.Bin`
            xlabel = "bin index $j$"
            ylabel = "bin index $i$"
            i = np.linspace(0, mg.N, mg.N)
            x, y = i, i
        z = K

        super().__init__(
            x, y, z, 
            xlabel=xlabel, ylabel=ylabel,
            *args, **kwargs
        )
