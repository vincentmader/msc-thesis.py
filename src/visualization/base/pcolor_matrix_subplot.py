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
        title: Optional[str] = None,
        xlabel: Optional[str] = None,
        ylabel: Optional[str] = None,
        symmetrized: bool = False,
        scales: tuple[str, str, str] = ("log", "log", "log"), # TODO rename?
        cmap: str = "Reds",
        z_limits: Optional[tuple[float, float]] = None,
        grid: Optional[bool] = False,
        xticks: Optional[Any] = None,
        yticks: Optional[Any] = None,
    ):
        for scale in scales:
            assert scale in ["lin", "log"], "Invalid scale."

        if len(z.shape) < 3:
            assert k is None, "2D matrix has no third index k."
            assert z.shape[0] == z.shape[1], "Non-cubic kernel shape."
        else:
            assert z.shape[0] == z.shape[1] == z.shape[2], "Non-cubic kernel shape."

        super().__init__(
            title=title,
            xlabel=xlabel,
            ylabel=ylabel,
        )
        self.x, self.y, self.z = x, y, z

        if symmetrized:
            if len(z.shape) == 2:
                self.z = 0.5 * (z + z.T) 
            elif len(z.shape) == 3:
                self.z = np.array([0.5 * (z + z.T) for z in z])
            else:
                raise Exception("")

        if len(z.shape) < 3:
            self.k = None
        else:
            self.k = k if k is not None else z.shape[2] // 2
        self.scales = scales
        self.cmap = cmap
        self.z_limits = z_limits
        self.grid = grid
        self.xticks = xticks
        self.yticks = yticks

    def draw(self, axes):
        k, x, y = self.k, self.x, self.y
        z = self.z if k is None else self.z[k]

        if self.scales[2] == "lin":
            if self.z_limits is None:
                vmin, vmax = z.min(), z.max()
            else:
                vmin = self.z_limits[0]
                vmax = self.z_limits[1]
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

        ax = axes[1]
        plt.sca(ax)
        plt.cla()
        self.im = plt.pcolormesh(
            x, y, z,
            cmap=self.cmap,
            norm=norm,
            rasterized=True
        )
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
        raise Exception("Not implemented.")
