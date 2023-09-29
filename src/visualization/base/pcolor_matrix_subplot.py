from typing import Optional, Any

import matplotlib.pyplot as plt
from matplotlib import colors
import numpy as np

from visualization.base import GridspecSubplot


class PcolorMatrixSubplot(GridspecSubplot):
    __slots__ = [
        "x", "y", "z", "k", "im", 
        "scales", "cmap", "z_limits", "grid",
        "xticks", "yticks",
    ]

    def __init__(
        self, 
        x: np.ndarray,
        y: np.ndarray,
        z: np.ndarray,
        k: Optional[int] = None,
        symmetrized: bool = False,
        scales: tuple[str, str, str] = ("log", "log", "log"), # TODO rename?
        cmap: str = "Reds",
        z_limits: Optional[tuple[float, float]] = None,
        grid: Optional[bool] = False,
        xticks: Optional[Any] = None,
        yticks: Optional[Any] = None,
        *args, **kwargs,
    ):
        # Assert valid scale specifications, & valid (cubic) matrix shape.
        for scale in scales:
            assert scale in ["lin", "log"], "Invalid scale."
        if len(z.shape) < 3:
            assert k is None, "2D matrix has no third index k, no need to specify it."
            assert z.shape[0] == z.shape[1], "Matrix shape should be cubic."
        else:
            assert z.shape[0] == z.shape[1] == z.shape[2], "Matrix shape should be cubic."

        # Initialize class instance.
        super().__init__(*args, **kwargs)
        self.scales, self.z_limits = scales, z_limits
        self.cmap, self.grid = cmap, grid
        self.xticks, self.yticks = xticks, yticks
        self.x, self.y, self.z = x, y, z

        # Symmetrize 3D matrix "layers".
        if symmetrized and len(z.shape) == 2:
            self.z = (z + z.T) / 2
        if symmetrized and len(z.shape) == 3:
            self.z = np.array([(z + z.T)/2 for z in z])

        # Define initial matrix "layer index".
        if len(z.shape) < 3:
            self.k = None
        else:
            self.k = k if k is not None else z.shape[2] // 2

    def draw(self, axes):
        k, x, y = self.k, self.x, self.y
        z = self.z if k is None else self.z[k]

        # Define norm.
        norm = self._pcolormesh_norm(z)

        ax = axes[1]
        plt.sca(ax)
        plt.cla()
        self.im = plt.pcolormesh(x, y, z, cmap=self.cmap, norm=norm, rasterized=True)
        plt.axis("scaled")
        plt.xlabel(self.xlabel)
        plt.ylabel(self.ylabel)
        plt.grid(self.grid)
        scale_x = self.scales[0]
        scale_y = self.scales[1]
        scale_x = "linear" if scale_x == "lin" else scale_x
        scale_y = "linear" if scale_y == "lin" else scale_y
        ax.set_xscale(scale_x)
        ax.set_yscale(scale_y)
        ax.format_coord = self.format_coord

        ax = axes[0]
        plt.sca(ax)
        plt.colorbar(self.im, cax=ax, orientation="horizontal", extend="both")
        plt.title(self.title)

    def update(self, k):
        self.k = k
        z = self.z[k]
        self.im.set_array(z.ravel())

    def format_coord(self, x, y):
        # TODO: Redefine.
        return f"{x}, {y}"

    def _pcolormesh_norm(self, z):
        if self.scales[2] == "lin":
            if self.z_limits is None:
                vmin, vmax = z.min(), z.max()
            else:
                vmin, vmax = self.z_limits
            smin, smax = np.sign(vmin), np.sign(vmax)
            if smin == smax or 0 in [smin, smax]:
                norm = colors.Normalize(vmin=vmin, vmax=vmax)
            else:
                norm = colors.TwoSlopeNorm(vmin=vmin, vcenter=0, vmax=vmax)
        elif self.scales[2] == "log":
            if self.z_limits is None:
                if z.min() == z.max() == 0:
                    vmin = 1e-20
                    vmax = 1e-7
                    print(f"Reverting to [{vmin}, {vmax}]")
                else:
                    vmin = z.min() if z.min() != 0 else z[z != 0].min()
                    vmax = z.max() if z.max() != 0 else z[z != 0].max()
            else:
                vmin = self.z_limits[0]
                vmax = self.z_limits[1]
            norm = colors.LogNorm(vmin=vmin, vmax=vmax)
        else:
            raise Exception("")
        return norm
