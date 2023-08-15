from typing import Any, Optional

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import LogLocator, LogFormatterMathtext

from kernel import Kernel
from dust import particle_radius_from_mass

# Features:
# 1.                                                                -> done
#   a) Plot a single 2D matrix.
#   b) Plot multiple 2D matrices. (1-3)
# 2.                                                                -> done
#   a) Plot with colorbar.
#   b) Plot without colorbar.
# 3.
#   a) Plot interactively with slider.
#   b) Plot non-interactively without slider.
# 4.
#   a) Plot kernel K_kij
#   b) Plot rel. vel. v_ij
#   c) Plot coll. rate C_ij
#   c) Plot coll. outcome prob. P_ij
# 5. 
#   a) Plot with lin. scaling in the pcolor plot & colorbar.
#   a) Plot with log. scaling in the pcolor plot & colorbar.
# 6. 
#   a) Plot with lin. scaling on the x- and y- axes.
#   a) Plot with log. scaling on the x- and y- axes.


class GridspecSubplot():
    __slots__ = ["data", "title", "xlabel", "ylabel"]

    def __init__(
        self, 
        data:       Any,
        title:      Optional[str]   = None,
        xlabel:     Optional[str]   = None,
        ylabel:     Optional[str]   = None,
    ):
        self.data = data
        self.title = title
        self.xlabel = xlabel
        self.ylabel = ylabel

    def draw(self, _):
        raise Exception("Not implemented.")


class KernelSubplot(GridspecSubplot):
    __slots__ = ["k"]

    def __init__(
        self, 
        kernel:     Kernel,
        k:          Optional[int]   = None,
        title:      Optional[str]   = None,
        xlabel:     Optional[str]   = None,
        ylabel:     Optional[str]   = None,
    ):
        super().__init__(
            kernel, 
            title=title,
            xlabel=xlabel,
            ylabel=ylabel,
        )
        self.k = k if k is not None else 25 #  TODO

    def draw(self, axes):
        kernel = self.data
        K = kernel.K_gain
        k = self.k
        K_k = K[k]

        mg = kernel.mg
        mc = mg.grid_cell_centers
        rho_s = 1600                          # TODO
        ac = particle_radius_from_mass(mc, rho_s)

        ax = axes[1]
        plt.sca(ax)
        im = plt.pcolor(
            ac, ac, 
            K_k, 
            cmap="Reds",
            norm=mpl.colors.LogNorm(
                vmin=1e-30, # vmin=data.min(),  TODO
                vmax=1e+15, # vmax=data.max(),  TODO
            ),
        )
        plt.axis("scaled")
        plt.xlabel(self.xlabel)
        plt.ylabel(self.ylabel)
        ax.set_xscale('log')
        ax.set_yscale('log')

        ax = axes[0]
        plt.sca(ax)
        plt.colorbar(im, cax=ax, orientation="horizontal")
        plt.title(self.title)


class BasePlot:
    __slots__ = ["fig"]

    def __init__(
        self,
        figsize: Optional[tuple[int, int]] = None,
    ):
        self.fig = plt.figure(figsize=figsize)

    def render(
        self,
        show_plot=True,
        save_plot=False,
        path_to_outfile=None,
    ):
        self.draw()
        if show_plot:
            plt.show()
        if save_plot:
            assert path_to_outfile is not None, "No savefile path specified."
            plt.savefig(path_to_outfile)
        plt.close()

    def draw(self):
        raise Exception("Not implemented.")


class GridspecPlot(BasePlot):
    __slots__ = ["gs", "axes", "subplots"]

    def __init__(
        self,
        subplots: list[GridspecSubplot],
    ):
        if len(subplots) == 1:
            figsize = (5, 6)
        elif len(subplots) == 2:
            figsize = (10, 6)
        elif len(subplots) == 3:
            figsize = (15, 6)
        else:
            raise Exception("")
        super().__init__(figsize=figsize)

        self.subplots = subplots

        gs_row_sizes = [1, 20]
        gs_col_sizes = [1] * len(subplots)
        gs_dimensions = (2, len(subplots))
        gs_locations = []
        for idx, subplot in enumerate(subplots):
            gs_location = (0, 1, idx, idx+1)
            gs_locations.append(gs_location)
            gs_location = (1, 2, idx, idx+1)
            gs_locations.append(gs_location)

        self.gs = self.fig.add_gridspec(
            nrows=gs_dimensions[0], 
            ncols=gs_dimensions[1],
            width_ratios=gs_col_sizes,
            height_ratios=gs_row_sizes,
            hspace=0.0,
            wspace=0.2,
        )

        self.axes = []
        for gs_location in gs_locations:
            row_i, row_f, col_i, col_f = gs_location
            foo = self.gs[row_i:row_f, col_i:col_f]
            ax = self.fig.add_subplot(foo)
            self.axes.append(ax)

    def draw(self):
        for idx, subplot in enumerate(self.subplots):
            axes = self.axes[2*idx: 2*idx+2]
            subplot.draw(axes)
