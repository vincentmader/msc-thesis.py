from typing import Optional

import matplotlib.pyplot as plt
from matplotlib import colors
import numpy as np

from .gridspec_subplot import GridspecSubplot


class PcolorMatrixSubplot(GridspecSubplot):
    __slots__ = ["x", "y", "z", "k", "im", "scales", "cmap"]

    def __init__(
        self, 
        x: np.ndarray,
        y: np.ndarray,
        z: np.ndarray,
        k: Optional[int] = None,
        title: Optional[str] = None,
        xlabel: Optional[str] = None,
        ylabel: Optional[str] = None,
        symmetrize: bool = False,
        scales: tuple[str, str, str] = ("log", "log", "log"), # TODO rename?
        cmap: str = "Reds",
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

        if symmetrize:
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

    def draw(self, axes):
        k, x, y = self.k, self.x, self.y
        z = self.z if k is None else self.z[k]

        if self.scales[2] == "lin":
            vmin, vmax = z.min(), z.max()
            smin, smax = np.sign(vmin), np.sign(vmax)
            if smin == smax or 0 in [smin, smax]:
                norm = colors.Normalize(vmin=vmin, vmax=vmax)
            else:
                norm = colors.TwoSlopeNorm(vmin=vmin, vcenter=0, vmax=vmax)
                self.cmap = "bwr"
        elif self.scales[2] == "log":
            vmin = z.min() if z.min() != 0 else z[z != 0].min()
            vmax = z.max() if z.max() != 0 else z[z != 0].max()
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
        scale_x = self.scales[0]
        scale_y = self.scales[1]
        scale_x = "linear" if scale_x == "lin" else scale_x
        scale_y = "linear" if scale_y == "lin" else scale_y
        ax.set_xscale(scale_x)
        ax.set_yscale(scale_y)

        ax = axes[0]
        plt.sca(ax)
        plt.colorbar(self.im, cax=ax, orientation="horizontal")
        plt.title(self.title)

    def update(self, k):
        self.k = k
        z = self.z[k]
        self.im.set_array(z.ravel())
