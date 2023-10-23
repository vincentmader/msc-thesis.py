from typing import Optional, Any

import matplotlib.pyplot as plt
from matplotlib import colors
import numpy as np

from visualization.base import GridspecSubplot

VMIN, VMAX = 1e-20, 1e-7


class PcolorMatrixSubplot(GridspecSubplot):
    __slots__ = [
        "x", "y", "z", "k", "im", 
        "axis_scales", "cmap", "z_limits", "grid",
        "xticks", "yticks",
    ]

    def __init__(
        self, 
        x:              np.ndarray,
        y:              np.ndarray,
        z:              np.ndarray,
        xticks:         Optional[Any]                   = None,
        yticks:         Optional[Any]                   = None,
        z_limits:       Optional[tuple[float, float]]   = None,
        k:              Optional[int]                   = None,
        grid:           bool                            = False,
        symmetrized:    bool                            = False,
        cmap:           str                             = "Reds",
        axis_scales:    tuple[str, str, str]            = ("log", "log", "log"),
        *args, **kwargs,
    ):
        # Assert valid scale specifications, & valid (cubic) matrix shape.
        for scale in axis_scales:
            assert scale in ["lin", "log"], "Invalid scale."
        if len(z.shape) < 3:
            assert k is None, "2D matrix has no third index k, no need to specify it."
            # assert z.shape[0] == z.shape[1], "Matrix shape should be cubic."
            assert z.shape[0] == z.shape[1], "Matrix shape should be square."
        else:
            assert z.shape[1] == z.shape[2], "Matrix layer shape should be square."

        # Initialize class instance.
        super().__init__(*args, **kwargs)
        self.x, self.y, self.z = x, y, z
        self.cmap, self.grid = cmap, grid
        self.xticks, self.yticks = xticks, yticks
        self.axis_scales, self.z_limits = axis_scales, z_limits

        # Symmetrize 3D matrix "layers".
        if symmetrized and (len(z.shape) == 2):
            self.z = (z + z.T) / 2
        if symmetrized and (len(z.shape) == 3):
            self.z = np.array([(z + z.T)/2 for z in z])

        # Define initial matrix "layer index".
        if len(z.shape) < 3:
            self.k = None
        else:
            self.k = k if k is not None else z.shape[2] // 2

    def draw(self, axes):
        k, x, y = self.k, self.x, self.y
        z = self.z if k is None else self.z[k]

        ax = axes[1]
        plt.sca(ax)
        norm = self._pcolormesh_norm(z)
        self.im = plt.pcolormesh(x, y, z, cmap=self.cmap, norm=norm, rasterized=True)
        plt.axis("scaled")
        scale_x = self.axis_scales[0]
        scale_y = self.axis_scales[1]
        ax.set_xscale("linear" if scale_x == "lin" else scale_x)
        ax.set_yscale("linear" if scale_y == "lin" else scale_y)
        ax.set_xlabel(self.xlabel)
        ax.set_ylabel(self.ylabel)
        ax.format_coord = self.format_coord
        ax.grid(self.grid)

        ax = axes[0]
        plt.sca(ax)
        plt.colorbar(self.im, cax=ax, orientation="horizontal", extend="both")
        plt.title(self.title)

    def update(self, k):
        self.k, z = k, self.z[k]
        self.im.set_array(z.ravel())

    def format_coord(self, x, y):
        return f"{x}, {y}"  # NOTE: Redefine this in child classes.

    def _pcolormesh_norm(self, z):
        if self.axis_scales[2] == "lin":
            vmin, vmax = (z.min(), z.max()) if self.z_limits is None else self.z_limits
            smin, smax = np.sign(vmin), np.sign(vmax)
            if smin == smax or 0 in [smin, smax]:
                return colors.Normalize(vmin=vmin, vmax=vmax)
            else:
                return colors.TwoSlopeNorm(vmin=vmin, vcenter=0, vmax=vmax)
        if self.axis_scales[2] == "log":
            if self.z_limits is None:
                if z.min() == z.max() == 0:
                    vmin, vmax = VMIN, VMAX
                    print(f"Reverting to [{vmin}, {vmax}]")
                else:
                    vmin = z.min() if z.min() != 0 else z[z != 0].min()
                    vmax = z.max() if z.max() != 0 else z[z != 0].max()
            else:
                vmin, vmax = self.z_limits[0], self.z_limits[1]
            return colors.LogNorm(vmin=vmin, vmax=vmax)
