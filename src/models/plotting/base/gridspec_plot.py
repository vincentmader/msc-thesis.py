from typing import Optional

from matplotlib.widgets import Slider

from .base_plot import BasePlot
from .gridspec_subplot import GridspecSubplot


class GridspecPlot(BasePlot):
    __slots__ = ["gs", "axes", "subplots", "slider", "slider_label"]

    def __init__(
        self,
        subplots: list[GridspecSubplot],
        add_slider: bool = False,
        slider_label: str = "k",
        figsize: Optional[tuple[int, int]] = None,
        gridspec_dimensions: Optional[tuple[int, int]] = None,
        *args, **kwargs
    ):
        figsize = self._choose_figsize(len(subplots)) if figsize is None else figsize
        super().__init__(figsize=figsize, *args, **kwargs)
        self.subplots = subplots
        self.slider_label = slider_label
        self._setup_axes(subplots, add_slider, gridspec_dimensions)

    def draw(self):
        self.fig.canvas.draw_idle()
        for idx, subplot in enumerate(self.subplots):
            axes = self.axes[2*idx: 2*idx+2]
            subplot.draw(axes)

    def _choose_figsize(
        self, 
        nr_of_subplots: int,
    ):
        figsizes = {
            1: (5, 6),
            2: (10, 6),
            3: (15, 6),
            4: (20, 6),
        }
        return figsizes[nr_of_subplots]

    def _setup_axes(
        self, 
        subplots: list[GridspecSubplot], 
        add_slider: bool,
        gridspec_dimensions: Optional[bool],
    ):
        hspace=0.0
        wspace=0.3

        if gridspec_dimensions == (2, 2):
            gs_dimensions = (4, 2)
            gs_row_sizes = [1, 20, 1, 20]
            gs_col_sizes = [20, 20]
            gs_col_sizes += [1] if add_slider else []
            gs_locations = []
            for idx, subplot in enumerate(subplots):
                gs_location = 2*(idx//2)+0, 2*(idx//2)+1, idx%2, idx%2+1
                # gs_location = ((2*idx)//2 + 0, idx//2 + 1, (2*idx)%2, (2*idx)%2+1)
                gs_locations.append(gs_location)
                gs_location = 2*(idx//2)+1, 2*(idx//2)+2, idx%2, idx%2+1
                # gs_location = ((2*idx)//2 + 1, idx//2 + 2, (2*idx)%2, (2*idx)%2+1)
                gs_locations.append(gs_location)
            hspace=0.25
            wspace=0.3

        else:
            nr_of_rows =  2
            nr_of_cols = len(subplots) + 1 if add_slider else len(subplots)
            gs_dimensions = (nr_of_rows, nr_of_cols)
            gs_row_sizes = [1, 20]
            gs_col_sizes = [20] * len(subplots)
            gs_col_sizes += [1] if add_slider else []
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
            hspace=hspace,
            wspace=wspace,
        )

        self.axes = []
        for gs_location in gs_locations:
            row_i, row_f, col_i, col_f = gs_location
            foo = self.gs[row_i:row_f, col_i:col_f]
            ax = self.fig.add_subplot(foo)
            self.axes.append(ax)

        if add_slider:
            self._setup_slider()

    def _setup_slider(self):
        k_max = max([s.z.shape[0] - 1 for s in self.subplots]) # TODO
        k_init = k_max // 2 + 1

        ax = self.fig.add_subplot(self.gs[-1])
        self.slider = Slider(
            ax=ax, 
            valmin=0, 
            valmax=k_max, 
            valstep=1, 
            valinit=k_init,
            label=self.slider_label, 
            orientation="vertical",
        )  # TODO
        self.slider.on_changed(self.update)
        self.axes.append(ax)

    def update(
        self, 
        k: int,
    ):
        for subplot in self.subplots:
            subplot.k = k
            subplot.update(k)
